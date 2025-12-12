module FancyData
__precompile__()

function help()
    println("""
    FancyData.jl — Utilities for Measurements, DataFrames, XML-based fit extraction
    ------------------------------------------------------------------------------

    Core measurement utilities:
      • val(x)                 – extract value() broadcasted
      • unc(x)                 – extract uncertainty() broadcasted
      • mes(x)                 – return [value, uncertainty]
      • mes2(x)                – formatted value(unc) string
      • mean2(x)               – mean with standard deviation
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

    """)
    return nothing
end

using Measurements, Printf, DataFrames, DelimitedFiles, LsqFit, Statistics, XMLDict

bin(arr,z) = [sum(arr[i:i+z-1]) for i in 1:z:length(arr)] # rebin spectrum, array-length has to be divisible by 'z'

val(x) = Measurements.value.(x)
unc(x) = Measurements.uncertainty.(x)
mes(x) = [val(x),unc(x)]

function mes2(x)
    try
		typeof(x) == Missing ? (return "") : nothing
        y,dy = mes(x)
        dy == 0.0 ? (return "$(round(y,sigdigits=6))") : nothing # exact value case

        i = -Int(floor(log10(dy))) # precision digit
        first_digits = round(dy*10.0^i,sigdigits=2)
        if first_digits < 3.0
            i+=1
        end
        # dy = round(Int, dy * 10^maximum([i,1]))
        dy = round(Int, dy * 10.0^i)
        y = round(y, digits=i)

        if isinteger(y) && dy % 10 == 0 # trailing zero case
            i -= 1
            dy ÷= 10
        end

        # if i < 0 #dy >= 30 # scientific notation
        #     c = 0
        #     while y >= 10.0
        #         y /= 10
        #         c += 1
        #     end
        #     y =  round(y, digits=c+i)
        #     return "$y($dy)e$c"
        # end

        return @sprintf("%.*f(%d)", i, y, dy)
    catch
        return x
    end
end

function wmean2(x, flag::Symbol=:none) # flag=:inex for internal and external errors
    x, dx = mes(x)
    n = length(x)
    n == 1 ? (return x[1], dx[1], dx[1]) : nothing
    w = [1/(dx[i])^2 for i in 1:n]
    xbar = sum([w[i] * x[i] for i in 1:n]) / sum([w[i] for i in 1:n])
    dx_extern = sqrt(1/(n-1) * sum([w[i] * (x[i] - xbar)^2 for i in 1:n]) / sum([w[i] for i in 1:n]))
    dx_intern = 1 / sqrt(sum([w[i] for i in 1:n]))
    if flag === :inex
        return xbar, dx_extern, dx_intern
    else
        return  measurement(wmean2(x)[1], maximum(wmean2(x)[2:3]))
    end
end

function mean2(x) # normaler mittelwert aber mit standardabweichung
    n = length(x)
    xbar = sum([val.(x[i]) for i in 1:n])/n
    dx = sqrt(sum([(val.(x[i])-xbar)^2 for i in 1:n] / (n * (n - 1))))
    # dx = maximum([dx, unc(mean(x))])
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