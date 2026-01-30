"""
# FancyData.jl

## Main Functions

- `mes()`: Converts a Measurement to a String with format `value(uncertainty)`
- `mean_std()`: Calculates the mean with standard deviation of a given array
- `wmean()`: Calculates the weighted mean of a given array of Measurements
- `tableDF()`: Prints a DataFrame to a TeX-formatted table
- `readDF()`, `writeDF()`: Reads a delimited file to a DataFrame and vice versa
- `readfits()`: Fitlistreader for HDTV-based XML-files
"""
module FancyData
__precompile__()

using Measurements, Printf, DataFrames, DelimitedFiles, LsqFit, Statistics, EzXML

val(x) = Measurements.value(x)
unc(x) = Measurements.uncertainty(x)

"""
    mes(x::Measurement)

Format a measurement as a string value(uncertainty).
Returns "" for `missing` and a plain number for exact values.
"""
mes(x) = x
mes(x::Missing) = ""
function mes(x::Measurement)

    # Ensure a Measurement was input
    try
        y,dy = val(x),unc(x)
    catch err 
        @warn "mes: failed to interpret input, returning original value" input=x exception=err
        return x
    end

    if dy == 0.0
        out = string(round(y, sigdigits=6))
        out[end-1:end] == ".0" && (out = out[1:end-2])
        return out 
    end

    scientific_notation = dy >= 30 && y != 0
    if scientific_notation
        exp = floor(Int, log10(abs(y)))
        y /= 10.0^exp
        dy /= 10.0^exp
    end

    precision_number = -floor(Int, log10(dy))
    first_digits = round(dy * 10.0^precision_number, sigdigits=2) # 10.0 is important don't change it to 10!
    first_digits < 3 && (precision_number += 1)

    y = round(y, digits=precision_number)
    dy = round(Int, dy * 10.0^precision_number)

    if isinteger(y) && dy % 10 == 0
        precision_number -= 1
        dy ÷= 10
    end

    scientific_notation && return @sprintf("%.*f(%d)e%d", precision_number, y, dy, exp)
    return @sprintf("%.*f(%d)", precision_number, y, dy)
end

"""
    wmean(X::AbstractVector{<:Measurement}, mode::Symbol=:max)

Calculates the weighted mean of `X`. 

Available modes:
- `:max`: (default): a `Measurement` with the larger of internal/external error
- `:both`: Returns external and internal uncertainties
"""
function wmean(X::AbstractVector{<:Measurement}; mode::Symbol=:max)
    mode in [:max,:both] || throw(ArgumentError("mode must be :max or :both (got :$mode)"))
    x,dx = val.(X),unc.(X)
    l = length(x)
    weights = 1 ./ (dx.^2) 
    xbar = sum(weights .* x) / sum(weights)
    dx_extern = sqrt(1 / (l-1) * sum(weights .* (x.-xbar).^2) / sum(weights))
    dx_intern = 1 / sqrt(sum(weights))
    l == 1 && (dx_extern = dx_intern)
    mode === :both && return (mean=xbar, ext=dx_extern, int=dx_intern)
    return  measurement(xbar, max(dx_extern,dx_intern))
end

"""
    mean_std(X::AbstractVector)

Returns Measurement of mean with standard deviation as uncertainty
"""
function mean_std(X::AbstractVector)
    any(x -> x isa Measurement, X) && @warn "input contains Measurements; uncertainties are ignored"
    X = val.(X)
    l = length(X)
    xbar = sum(X) / l
    dx = sqrt(sum((X.-xbar).^2 / (l * (l - 1))))
    return measurement(xbar, dx)
end

function parse_measurement(s::AbstractString)
    # Regex with optional minus, optional decimal, optional parentheses, optional exponent
    m = match(r"^(-)?(\d+)(?:\.(\d+))?(?:\((\d+)\))?(?:e([+-]?\d+))?$", s)

    if m === nothing
        # Not a valid measurement string -> return zeros and no sign
        return (sign = "", int_digits = 0, frac_digits = 0, par_digits = 0, exp_digits = 0)
    end

    sign, intpart, fracpart, parpart, expart = m.captures

    return (
        sign = sign === nothing ? "" : "-",
        int_digits  = intpart  === nothing ? 0 : length(intpart),
        frac_digits = fracpart === nothing ? 0 : length(fracpart),
        par_digits  = parpart  === nothing ? 0 : length(parpart),
        exp_digits  = expart   === nothing ? 0 : length(expart)
    )
end

"""
    tableDF(DF::DataFrame)

Returns LaTeX table from a DataFrame.
"""
function tableDF(DF)
    l,w = size(DF)
    DF = mes.(DF)
    lengths = [maximum(length.([(typeof(DF[j,i])==Missing) ? 0 : string(DF[j,i]) for j in 1:l])) for i in 1:w]

    ## SIUnitx columns
    bunchastrings=["\\begin{table}[]\n\t\\centering\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\t\\begin{tabular}{\n"]
    for j in 1:w 
        a0 = any(parse_measurement("$(DF[i,j])").sign == "-" for i in 1:l) ? "-" : ""
        a1 = maximum(parse_measurement("$(DF[i,j])").int_digits  for i in 1:l)
        a2 = maximum(parse_measurement("$(DF[i,j])").frac_digits for i in 1:l)
        a3 = maximum(parse_measurement("$(DF[i,j])").par_digits  for i in 1:l)
        a4 = maximum(parse_measurement("$(DF[i,j])").exp_digits  for i in 1:l)
        a = "S[table-format=$a0$a1.$a2($a3)e$a4]"
        if a == "S[table-format=0.0(0)e0]"
            push!(bunchastrings,"\t\tl\n")
        else
            push!(bunchastrings,"\t\t$a\n")
        end
    end
    push!(bunchastrings,"\t}\n")

    push!(bunchastrings,"\t\t{"*(join(names(DF),"}  &  {")*"} \\\\ \\midrule\n")) #Kopfzeile
    for j in 1:l
        vals = [rpad("&  $(DF[j,i])",lengths[i]+5) for i in 1:w]
        vals[1] = rpad("\t\t$(DF[j,1])",lengths[1]+2)
        push!(vals,"\\\\\n")
        push!(bunchastrings,join(vals))
    end
    push!(bunchastrings,"\t\\end{tabular}\n\\end{table}")
    return join(bunchastrings)
end 

"""
    writeDF(file_out::String, DF::DataFrame; delim::Char='\\t')

Saves a DataFrame as a delimted file with default delimiter `\\t`.
"""
writeDF(out,DF;delim='\t') = writedlm(out,Iterators.flatten(([names(DF)],eachrow(DF))), delim)

"""
    readDF(file_in::String; delim::Char='\\t')

Reads a delimited file to a DataFrame with default delimiter `\\t` and attempts basic parseing of Measurements and Missing fields.
"""
function readDF(path_to_DF; delim='\t')
    data = readdlm(path_to_DF, delim)
    names = Symbol.(data[1,:])
    matrix = data[2:end,:]

    function transform_element(val)
        if val isa AbstractString
            if val == "missing"
                val = missing
            else
                m = tryparse(Measurement{Float64}, val)
                return isnothing(m) ? val : m
            end
        else
            return val
        end
    end
    return DataFrame(transform_element.(matrix), names)
end 

"""
    readfits(XML_file::String; mode::Symbol=:peak)

Reads an HDTV file and returns DataFrame

Available modes:
- `:peak`: position, width and volume of fitted peaks
- `:sub`: position, width and volume of background-subtracted integrals
- `:param`: peak, region and background marker of fits (if available)

All returned values are calibrated if possible with uncalibrated values as a fallback.

If uncertainties are present, the values are given as measurements.
"""
function readfits(file::String; mode::Symbol=:peak)
    document = readxml(file)
    fits = findall("//fit", document)

    # Case 1: mode is :peak or :sub => pos, wid and vol
    if mode in [:peak, :sub]
        DF = DataFrame(
            ID    = Float64[],
            pos   = Measurement{Float64}[],
            width = Measurement{Float64}[],
            vol   = Measurement{Float64}[]
        )
    
        for (i,fit) in enumerate(fits)
            
            # Select Type 
            if mode == :peak
                peak_nodes = findall("peak", fit)
            else
                peak_nodes = findall("integral[@integraltype='sub']", fit)
            end
    
            for (j,peak_node) in enumerate(peak_nodes)
                next_row = DataFrame(
                    ID    = [0.0],
                    pos   = [measurement(0.0, 0.0)], 
                    width = [measurement(0.0, 0.0)], 
                    vol   = [measurement(0.0, 0.0)]
                )
        
                # Use calibrated values if available, uncalibrated otherwise
                calibrated_node = findfirst("cal", peak_node)
                isnothing(calibrated_node) && (calibrated_node = findfirst("uncal", peak_node))

                for quantity in ["pos", "width", "vol"]
                    quantity_node = findfirst(quantity, calibrated_node)
                    value_node = findfirst("value", quantity_node)
                    error_node = findfirst("error", quantity_node)
                    value = !isnothing(value_node) ? parse(Float64, nodecontent(value_node)) : 0.0
                    error = !isnothing(error_node) ? parse(Float64, nodecontent(error_node)) : 0.0
                    next_row[1, quantity] = measurement(value, error)
                end
    
                next_row[1, "ID"] = i+(j-1)/10
                
                append!(DF, next_row)    
            end
            
        end

    # Case 2: mode is :param => parameter-markers
    elseif mode == :param
        DF = DataFrame(
            ID           = Float64[],
            peakMarker   = Vector{Float64}[],
            regionMarker = Tuple{Float64, Float64}[],
            bgMarker     = Vector{Tuple{Float64, Float64}}[]
        )

        for (i,fit) in enumerate(fits)

            next_row = DataFrame(
                ID           = [i],
                peakMarker   = [[]],
                regionMarker = (0.0, 0.0),
                bgMarker     = [Tuple{Float64, Float64}[]]
            )

            # Pear Markers
            peak_nodes = findall("peakMarker", fit)
            if !isempty(peak_nodes)
                for (j,peak_node) in enumerate(peak_nodes)
                    position_node = findfirst("position", peak_node)

                    calibrated_node = findfirst("cal", position_node)
                    isnothing(calibrated_node) && (calibrated_node = findfirst("uncal", position_node))
    
                    push!(next_row[1, "peakMarker"], parse(Float64, nodecontent(calibrated_node)))
                end
            end

            # Region Markers
            region_node = findfirst("regionMarker", fit)
            begin_node = findfirst("begin", region_node)
            end_node = findfirst("end", region_node)

            calibrated_begin_node = findfirst("cal", begin_node)
            calibrated_end_node = findfirst("cal", end_node)
            isnothing(calibrated_begin_node) && (calibrated_begin_node = findfirst("uncal", begin_node))
            isnothing(calibrated_end_node) && (calibrated_end_node = findfirst("uncal", end_node))

            region_begin = parse(Float64, nodecontent(calibrated_begin_node))
            region_end = parse(Float64, nodecontent(calibrated_end_node))

            next_row[1, "regionMarker"] = (region_begin, region_end)
            
            # Background Markers
            bg_nodes = findall("bgMarker", fit)
            if !isempty(bg_nodes)
                for (j,bg_node) in enumerate(bg_nodes)
                    begin_node = findfirst("begin", bg_node)
                    end_node = findfirst("end", bg_node)

                    calibrated_begin_node = findfirst("cal", begin_node)
                    calibrated_end_node = findfirst("cal", end_node)
                    isnothing(calibrated_begin_node) && (calibrated_begin_node = findfirst("uncal", begin_node))
                    isnothing(calibrated_end_node) && (calibrated_end_node = findfirst("uncal", end_node))

                    bg_begin = parse(Float64, nodecontent(calibrated_begin_node))
                    bg_end = parse(Float64, nodecontent(calibrated_end_node))
                    
                    push!(next_row[1, "bgMarker"], (bg_begin, bg_end))
                end
            end

            append!(DF, next_row)
        end
        
    else
        error("Mode must be :peak, :sub or :param")
    end
    
    return DF
end

export val, unc, mes, wmean, mean_std, tableDF, writeDF, readDF, readfits

end
