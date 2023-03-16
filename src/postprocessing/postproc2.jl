"""
    compute_intensity(z, I, κ, ϵ)
    simple integration alog z

    z - positions
    I(z) - Intensity
    κ(z) - absorption
    ϵ(z) - emission

    returns: intensity
"""
function compute_intensity(z, I, κ, ϵ)
    dz = 100.0
    exp_κ = exp(-κ * dz)
    a     = ϵ / κ * (1.0 - exp_κ)
    I_    = zeros(Float64, size(z,1))
    I_[1] = I[1]
    limit = 0.01
    ee   = 0.0
    for i in 1:size(z,1)-1
        ee = ee + ϵ[i] * dz
        if κ[i] * dz < limit
            I_[i+1] = I_[i] * exp_κ[i] + ϵ[i] * dz
        else
            I_[i+1] = I_[i] * exp_κ[i] + a[i]
        end
    end
    return I_
end


"""
# iN = 1, iθ = 1, NCO2 =  4.00000e+02, θ =  0.00000e+00

iz =   0, z =  0.00000e+00,  NCO2 =  4.00000e+02, θ =  0.00000e+00, I =  1.10178e+02, T =  2.88000e+02, N =  2.54761e+25, ϵ =  0.00000e+00, κ =  0.00000e+00, ΔλL =  0.00000e+00, ΔλD =  0.00000e+00
"""
function decode_line(line)
    l1 = split(line, ",")
    keys   = Vector{String}(undef,0)
    values = Vector{Float64}(undef,0)
    for (i, a) in enumerate(l1)
        b = split(a, "=")
        push!(keys, b[1])
        push!(values, parse(Float64, b[2]))
    end
    keys, values
end


"""[summary]

Arguments:
    id {[type]} -- [description]
    out_root {[type]} -- [description]
"""
function read_data(out_root)
    parameter = JSON.parsefile(joinpath(out_root, "input.json"))#; dicttype=Dict, inttype=Int64, use_mmap=true)

    lines = readlines(joinpath(out_root, "log.log"), keep=false) # dont keeep \n

    ls = split(lines[1])

    # 2  3  700  10
    # # iN = 1, iθ = 1, NCO2 =  4.00000e+02, θ =  0.00000e+00

    nN   = parse(Int64, ls[1])     # nb NCO2
    nθ   = parse(Int64, ls[2])     # nb theta
    nz   = parse(Int64, ls[3])+1   # nb z, there are two iz = 0 rows in log.log file
    ncol = parse(Int64, ls[4])     # nb columns
    #        0     1  2  3  4  5  6  7  8    9
    #    iz, z, NCO2, θ, I, T, N, ϵ, κ, ΔλL, ΔλD
    #iz =   0, z =  0.00000e+00,  NCO2 =  4.00000e+02, θ =  0.00000e+00, I =  1.10178e+02, T =  2.88000e+02, N =  2.54761e+25, ϵ =  0.00000e+00, κ =  0.00000e+00, ΔλL =  0.00000e+00, ΔλD =  0.00000e+00

    # 3  2  701 10
    M = zeros(Float64, nN, nθ, nz, ncol+1)
    iN = 1
    iθ = 1
    inext = 0
    for (i,line) in enumerate(lines[2:end])
        #println(line, "  ", line[1] == "#" )
        if line[1] == '#' # indicates next run
            k,v = decode_line(line[2:end])
            iN = floor(Int64, v[1])
            iθ = floor(Int64, v[2])
            inext = i+1
        else
            k,v = decode_line(line)
            iz = floor(Int64, v[1])
            M[iN, iθ, iz, :] = v[2:end]
            if i == inext
                @infoe @sprintf("%d  %d  %d  %s", iN, iθ, iz, v[2:end])
            end
        end
    end
    M, parameter
end
