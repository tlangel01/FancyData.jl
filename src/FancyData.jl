"""
FancyData.jl — Utilities for Measurements, DataFrames, XML-based fit extraction

Core measurement utilities:
    • val(x)                 – extract value() broadcasted
    • unc(x)                 – extract uncertainty() broadcasted
    • mes(x)                 – formatted value(unc) string
    • mean_std(x)            – mean with standard deviation
    • wmean(x; :inex)        – weighted mean as Measurement

DataFrame utilities:
    • tableDF(df)            – convert table to LaTeX SIunitx format
    • writeDF(path, df)      – write DataFrame with header
    • readDF(path; delim)    – read tabular file into DataFrame

Spectrum utilities:
    • bin(arr, z)            – rebin array into groups of z

Parenthesis utilities:
    • p2v(str)               – parse "1.23(4)" → measurement
    • p2c(in, out, cols)     – split column into value/error columns in new file

XML fit extractors:
    • readvals(xml, col; cal, type)
    • readfits(xml; cal, type) → DataFrame(:pos, :wid, :vol)

"""
module FancyData
__precompile__()

using Measurements, Printf, DataFrames, DelimitedFiles, LsqFit, Statistics, XMLDict

val(x::Measurement) = Measurements.value(x)
unc(x::Measurement) = Measurements.uncertainty(x)

"""
    mes(x::Measurement)

Format a `Measurement` as a string `value(uncertainty)`.
Returns `""` for `missing` and a plain number for exact values.
"""
function mes(x::Measurement)
    x === missing && return ""

    # Ensure a Measurement was input
    try
        y,dy = val(x),unc(x)
    catch err 
        @warn "mes: failed to interpret input, returning original value" input=x exception=err
        return x
    end

    dy == 0.0 && return string(round(y, sigdigits=6))

    # Determine precision
    i = -floor(Int, log10(dy))
    first_digits = round(dy * 10^i, sigdigits=2)
    first_digits < 3 && (i += 1)

    y = round(y, digits=i)
    dy = round(Int, dy * 10^i)

    # Handle trailing zero uncertainty
    if isinteger(y) && dy % 10 == 0
        i -= 1
        dy ÷= 10
    end

    return @sprintf("%.*f(%d)", i, y, dy)
end

"""
    wmean(X::AbstractVector{<:Measurement}, mode::Symbol=:max)

Calculates the weighted mean of `X`. 

Returns:
- `:max` (default): a `Measurement` with the larger of internal/external error
- `:both`: `(value, external_error, internal_error)`
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
    mode === :both && return xbar, dx_extern, dx_intern
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
        # Not a valid measurement string → return zeros and no sign
        return (sign = "", int_digits = 0, frac_digits = 0, par_digits = 0, exp_digits = 0)
    end

    sign, intpart, fracpart, parpart, expart = m.captures

    return (
        sign = sign === nothing ? "" : "-",
        int_digits = intpart === nothing ? 0 : length(intpart),
        frac_digits = fracpart === nothing ? 0 : length(fracpart),
        par_digits = parpart === nothing ? 0 : length(parpart),
        exp_digits = expart === nothing ? 0 : length(expart)
    )
end

function tableDF(DF)
    l,w = size(DF)
    DF = mes2.(DF)
    lengths = [maximum(length.([(typeof(DF[j,i])==Missing) ? 0 : string(DF[j,i]) for j in 1:l])) for i in 1:w]

    ## SIUnitx columns
    bunchastrings=["\\begin{table}[]\n\t\\centering\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\t\\begin{tabular}{\n"]
    for j in 1:w 
        ## Finde die Längen von vor/nachkommastellen und fehlern (sins against humanity)
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

writeDF(out,DF) = writedlm(out,Iterators.flatten(([names(DF)],eachrow(DF))),'\t')
readDF(path_to_DF;delim='\t') = DataFrame(readdlm(path_to_DF,delim)[2:end,:], strip.(readdlm(path_to_DF,delim)[1,:]))

function p2v(str) # parenthesis to value
    try
        return mes(measurement("$str"))
    catch
        return str 
    end
end

function p2c(infile,outfile,cols) # parentheses to column
    InArr = DataFrame(load(infile))
    Header = ["" for i in 1:(length(InArr[1,:]) + length(cols))]
    w_old = length(InArr[1,:])
    w = w_old + length(cols)
    l = length(InArr[:,1])
    s = 0
    OutArr = Array{Any}(undef,l,w)
    for c in 1:w_old
        c in cols ? (OutArr[:,c+s:(c+s+1)] = reduce(vcat,transpose.(p2v.(InArr[:,c]))); Header[c+s:c+s+1] = [names(InArr)[c],"Δ"*names(InArr)[c]]; s+=1) : (OutArr[:,c+s] = InArr[:,c]; Header[c+s]=names(InArr)[c])
    end
    save(outfile,DataFrame(OutArr,Header))
end

function readvals(file,col;cal=true,type="peak")
    A = []
    cal ?  cal = "cal" : cal = "uncal"
    for i in 1:length(file["fit"])
        if typeof(file["fit"][i][type]) == Vector{XMLDict.XMLDictElement}
            push!(A,[measurement.(parse.(Float64,file["fit"][i][type][j][cal][col]["value"]),parse.(Float64,file["fit"][i][type][j][cal][col]["error"])) for j in 1:length(file["fit"][i][type])])
        else
            push!(A,measurement.(parse.(Float64,file["fit"][i][type][cal][col]["value"]),parse.(Float64,file["fit"][i][type][cal][col]["error"])))
        end
    end
    return vcat(A...)
end

function readfits(file;cal=true,type="peak")
    str=parse_xml(read(file,String))
    pos=readvals(str,"pos", cal=cal, type=type)
    wid=readvals(str,"width", cal=cal, type=type)
    vol=readvals(str,"vol", cal=cal, type=type)
    return (DataFrame(([pos wid vol]),[:pos,:wid,:vol]))
end

export measurement, bin, val, unc, mes, mes2, wmean2, wmean, mean2, tableDF, writeDF, readDF, p2v, p2c, readvals, readfits




end # module 