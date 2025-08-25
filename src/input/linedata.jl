using PhysConst
using SimpleLog
using Printf
using StaticArrays
using Interpolations
using DataFrames
using CSV
using LoopVectorization
using Bumper


@inline function S_T(S21r, E1, E2, β, βr, QT, QTr)
    ΔE21 = E2-E1
    S21r * QTr/QT * exp(-E1*(β+βr)) * (1.0-exp(-ΔE21*β)) / (1.0-exp(-ΔE21*βr))
end

struct LineData
    species :: Symbol
    iso   :: Vector{Int64}
    λ210  :: Vector{Float64}
    ΔE21  :: Vector{Float64}
    E1    :: Vector{Float64}
    E2    :: Vector{Float64}
    A21   :: Vector{Float64}
    B21   :: Vector{Float64}
    B12   :: Vector{Float64}
    g2    :: Vector{Float64}
    g1    :: Vector{Float64}
    S21r  :: Vector{Float64} # at TREF
    γair  :: Vector{Float64} # at TREF
    γself :: Vector{Float64} # at TREF
    nair  :: Vector{Float64}
    δair  :: Vector{Float64} # at TREF
end

function get_nλl(ld::LineData)
    length(ld.A21)
end

function load_hitran_data(hitran_out, λmin, λmax, iso_max)
    # molec_id, local_iso_id, nu, sw, a, gamma_air, gamma_self, elower, n_air, delta_air, gp, gpp
    df0 = CSV.read(hitran_out, DataFrame)

    nu_m = df0[!,:nu] ./ _cm
    λ  = 1.0 ./ nu_m
    λ1,λ2 = extrema(λ)
    @infoe @sprintf("λ1 = %8.2e, λ2 = %8.2e, λmin = %8.2e, λmax = %8.2e", λ1, λ2, λmin, λmax)

    νmax = _cm/λmin
    νmin = _cm/λmax
    ids(x) = @. ( (x >= νmin) && (x < νmax) )
    df = df0[ids(df0.nu),:]

    iso_min = 1
    iso_max = min(maximum(df[!,:local_iso_id]), iso_max)
    isos(x) = @. ( (x >= iso_min) && (x <= iso_max) )
    df = df[isos(df.local_iso_id),:]

    sort!(df, :nu, rev=true)

    mid   = df[!,1]
    iso   = df[!,2]
    ν21   = df[!,3]                                                  # [1/cm]
    S21r  = df[!,4]                                                  # [cm^−1/(molecule ⋅ cm^−2)] = [cm / molecule], the spectral line intensity 
    A21   = df[!,5]                                                  # [[1/s]
    γair  = df[!,6]                                                  # [1 / (cm * atm)] = 1 / (1.0e-2 m * 1.01325e5 Pa)]
    γself = df[!,7]                                                  # [1/(cm * atm)]
    ν1    = df[!,8]                                                  # [1/cm]
    nair  = df[!,9]
    δair  = df[!,10]                                                 # [1/(cm * atm)]
    g2    = df[!,11]
    g1    = df[!,12]

    hc    = c_h*c_c

    ν1_m  = ν1  ./ _cm                                               # [1/m]
    ν21_m = ν21 ./ _cm                                               # [1/m]
    λ210  = 1.0 ./ ν21_m                                             # [m]
    E1    = hc  .* ν1_m                                              # [J]
    ΔE21  = hc  ./ λ210                                              # [J] 
    E2    = E1  .+ ΔE21                                              # [J]

    γair  = γair  ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal
    γself = γself ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal
    δair  = δair  ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal

    S21r = S21r .* _cm                                               # [m]

    # Einstein coefficient of induced emission
    B21 = @. A21 * λ210^3 / (8π * c_h)                               # [m^3 / s / Js] = [m^3 / J / s^2]
    #B21 = @. A21 * λ210^3 / (2.0 * c_h)                               # [m^3 / s / Js] = [m^3 / J / s^2]
    # Einstein coefficient of absorption
    B12 = @. g2 / g1 * B21;

    iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair
end

function hitran_to_hdf5(species, hitran_out, hdf5, hdf5c, λmin, λmax, iso_max)
    iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair = load_hitran_data(hitran_out, λmin, λmax, iso_max)

    datasets = Dict(string(species) => (
        ("iso"     , iso   ),
        ("λ210"    , λ210  ),
        ("ΔE21"    , ΔE21  ),
        ("E1"      , E1    ),
        ("E2"      , E2    ),
        ("A21"     , A21   ),
        ("B21"     , B21   ),
        ("B12"     , B12   ),
        ("g2"      , g2    ),
        ("g1"      , g1    ),
        ("S21r"    , S21r  ),
        ("γair"    , γair  ),
        ("γself"   , γself ),
        ("nair"    , nair  ),
        ("δair"    , δair  )
    ))

    @infoe hdf5, string(species)

    save_arrays_to_hdf5(hdf5, datasets)

    data = Matrix{Float64}(undef, 14, length(λ210))
    data[ 1,:] = λ210  
    data[ 2,:] = ΔE21  
    data[ 3,:] = E1    
    data[ 4,:] = E2    
    data[ 5,:] = A21   
    data[ 6,:] = B21   
    data[ 7,:] = B12   
    data[ 8,:] = g2    
    data[ 9,:] = g1    
    data[10,:] = S21r  
    data[11,:] = γair  
    data[12,:] = γself 
    data[13,:] = nair  
    data[14,:] = δair

    save_arrays_to_hdf5(hdf5c, Dict(string(species) => (("iso", iso), ("data", data))))
end

"""
    LineData(species, hitran_out, λmin, λmax, iso_max)    

    transform HITRAN data from hitran_file format (joinpath(data_dir, string(hitran_file, ".out")))
    to .npy (joinpath(input_dir, string(hitran_file, ".npy") and .hdf5 joinpath(input_dir, string(hitran_file, ".hdf5") formats, 
    using only values used in the computations

    
    species, hitran_out, λmin, λmax, iso_max = :CO2, CO2out, par[:λmin], par[:λmax], length(mdCO2.iso_a)
"""
function LineData(hdf5_path)
    iso_data = load_arrays_from_hdf5(hdf5_path)
    species = collect(keys(iso_data))[1]
    
    ld = if occursin("compact", hdf5_path)
        data  = iso_data[species]["data"]
        iso   = iso_data[species]["iso"]
        λ210  = data[ 1,:]
        ΔE21  = data[ 2,:]
        E1    = data[ 3,:]
        E2    = data[ 4,:]
        A21   = data[ 5,:]
        B21   = data[ 6,:]
        B12   = data[ 7,:]
        g2    = data[ 8,:]
        g1    = data[ 9,:]
        S21r  = data[10,:]
        γair  = data[11,:]
        γself = data[12,:]
        nair  = data[13,:]
        δair  = data[14,:]
        LineData(Symbol(species), iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair)
    else
        iso_data[species]["data"]
        iso   = iso_data["iso"]
        λ210  = iso_data["λ210"]                                              # m
        ΔE21  = iso_data["ΔE21"]                                              # J        
        E1    = iso_data["E1"]                                                # J
        E2    = iso_data["E2"]                                                # J
        A21   = iso_data["A21"]                                               # 1/s
        B21   = iso_data["B21"]                                               # m^3 / (J * s^2)
        B12   = iso_data["B12"]                                               # m^3 / (J * s^2)
        g2    = iso_data["g2"]                                                # 
        g1    = iso_data["g1"]                                                # 
        S21r  = iso_data["S21r"]                                              # m
        γair  = iso_data["γair"]                                              # 1 / (m * Pa)
        γself = iso_data["γself"]                                             # 1 / (m * Pa)
        nair  = iso_data["nair"]                                              # 
        δair  = iso_data["δair"]                                              # 1 / (m * Pa)
        LineData(Symbol(species), iso, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, S21r, γair, γself, nair, δair)
    end
    #@infoe @sprintf("nλl = %d", length(ld.λ210))
    ld
end

"""
    compute_line_emission_and_absorption_iλ(ld::LineData, Qref, Qiso, miso, c, T, N, p, iλ)
"""
@inline function compute_line_emission_and_absorption_iλ(line_data::LineData, Qref, Qiso, miso, aiso, cspech, T, N, p, iλl)
    dΩ = 1.0
    β  = 1.0/(c_kB * T)
    βr = 1.0/(c_kB * TREF)

    iso   = line_data.iso[iλl]                                               # 
    λ210  = line_data.λ210[iλl]                                              # m
    E1    = line_data.E1[iλl]                                                # J
    E2    = line_data.E2[iλl]                                                # J
    A21   = line_data.A21[iλl]                                               # 1/s
    B21   = line_data.B21[iλl]                                               # m^3 / (J * s^2)
    B12   = line_data.B12[iλl]                                               # m^3 / (J * s^2)
    g2    = line_data.g2[iλl]                                                # 
    g1    = line_data.g1[iλl]                                                # 
    S21r  = line_data.S21r[iλl]                                              # m
    γair  = line_data.γair[iλl]                                              # 1 / (m * Pa)
    γself = line_data.γself[iλl]                                             # 1 / (m * Pa)
    nair  = line_data.nair[iλl]                                              # 
    δair  = line_data.δair[iλl]                                              # 1 / (m * Pa)

    #      species concentartion at h, isotop abundance, gas density
    Niso = cspech * aiso[iso] * N                                            #  [1/m^3]

    # pressure shift
    λ21 = λ210 / (1.0 + δair * λ210 * p)

    # Lorentzian (pressure-broadened) HWHM, γ(p,T) 
    # γpT = (TREF/T)^nair * (γair * (p - pself) + γself*pself)
    γp = (TREF/T)^nair * (γair * p * (1.0 - cspech) + γself * p * cspech)      # [1/m]
    ΔλL = λ21^2 * γp                                                 # [m]

    # Doppler broadening
    ΔλG = sqrt(2.0 * c_kB * T / miso[iso]) / c_c * λ21

    # occupation numbers
    N1  = g1 * exp(- E1 * β) / Qiso[iso] * Niso
    N2  = g2 * exp(- E2 * β) / Qiso[iso] * Niso

    # emission [W/m^2]
    ϵ = hc/λ21 * N2 * A21 * dΩ / (4.0 * π)
    # ϵ = c_h * λ21 * A21 * dΩ/(4.0*π) *  N2                          # [J / (m^2 * sr)]
    # ϵ * f(λ) * dz                                                   # [J / (m^2 * sr) / m * m] 

    # absorption coefficient [1]
    κ = c_h * λ21 / c_c * (N1 * B12 - N2 * B21)                       # Js * m * s/m / m^3 * m^3/(J*s^2) = 1

    S21  = S_T(S21r, E1, E2, β, βr, Qiso[iso], Qref[iso]) * Niso      # [1/m^2]  
    #κ2   = S21 * λ210^2                                               # [1]

    SVector{13, Float64}(
        Float64(iso),
        miso[iso],
        aiso[iso],
        Qiso[iso],
        S21,
        λ21,
        γp,
        ΔλL,
        ΔλG,
        N1,
        N2,
        ϵ,
        κ)
end

@doc raw"""
    compute_lines_emission_and_absorption(moleculardata, linedata, Nmolecules, T, N, p)

    Compute emission and absorption coefficients of the lines

    T - temperature
    N - atmosphere density
    p - pressure
    NCO2 - CO2 concentration
No4
"""
function compute_lines_emission_and_absorption!(linedata::Vector{SVector{13, Float64}}, par, line_data::LineData, Qref, Qiso, miso, aiso, cspech, T, N, p)
    iλl = argmin(line_data.E1)
    iλl = 48512
    Threads.@threads for iλl in eachindex(line_data.λ210)
        linedata[iλl] = compute_line_emission_and_absorption_iλ(line_data, Qref, Qiso, miso, aiso, cspech, T, N, p, iλl)
    end
end

"""
    sum over all lines using their line shape

    linedata = linedata_dict[spec] 
"""
function sum_over_lines!(par, λb, linedata, ϵbt, κbt, ϵb, κb, int_f, int_ft, pr::Profile)

    f_Δλ_factor = par.r.f_Δλ_factor
    f_adapt     = par.r.f_adapt

    Dλ   = λb[end] - λb[1]
    Δλ   = λb[2]   - λb[1]
    nλb  = length(λb)

    # 1    2          3          4          5    6    7   8    9    10  11  12 13
    # iso, miso[iso], aiso[iso], Qiso[iso], S21, λ21, γp, ΔλL, ΔλG, N1, N2, ϵ, κ 

    λ21  = [linedata[i][6]  for i in eachindex(linedata)]
    ΔλLh = [linedata[i][8]  for i in eachindex(linedata)] .* 0.5
    ΔλGh = [linedata[i][9]  for i in eachindex(linedata)] .* 0.5
    ϵ    = [linedata[i][12] for i in eachindex(linedata)]
    κ    = [linedata[i][13] for i in eachindex(linedata)]
    nλl  = length(λ21)  

    nthreads = Threads.nthreads()
    for i in eachindex(ϵbt)
        fill!(ϵbt[i], 0.0)
        fill!(κbt[i], 0.0)
    end
    
    nthreads = Threads.nthreads()
    nf = floor(Int64, maximum((ΔλLh + ΔλGh)) * f_Δλ_factor / Δλ) * 2 + 10
    update_profile(pr, nf, nthreads, nλl)
    
    Threads.@threads for iλl in 1:nλl 
        tid = Threads.threadid()

        iλb = floor(Int64, (λ21[iλl] - λb[1]) / Dλ * Float64(nλb-1)) + 1
        δiλ = floor(Int64, (ΔλLh[iλl] + ΔλGh[iλl]) * f_Δλ_factor / Δλ)
        iλm = max(1, iλb - δiλ)
        iλp = min(nλb, iλb + δiλ + 1)

        λrange = @view λb[iλm:iλp]
        nλrange = length(λrange)
        if nf < nλrange
            @infoe iλl, nf, nλrange
        end
        
        ft = pr.ft[tid]
        voigt!(ft, λrange, λb[iλb], ΔλLh[iλl], ΔλGh[iλl], f_adapt)
        int_ft[tid][iλl] = sum(ft[1:nλrange])*Δλ;

        @turbo for iλ in iλm:iλp 
            jλ = iλ-iλm+1
            ϵbt[tid][iλ] = ϵbt[tid][iλ] + ft[jλ] * ϵ[iλl]
        end
        @turbo for iλ in iλm:iλp 
            jλ = iλ-iλm+1
            κbt[tid][iλ] = κbt[tid][iλ] + ft[jλ] * κ[iλl]
        end
    end

    fill!(ϵb, 0.0)
    fill!(κb, 0.0)
    for tid in 1:nthreads
        @. ϵb[:] += ϵbt[tid][:]
        @. κb[:] += κbt[tid][:]
    end

    fill!(int_f, 0.0)
    for tid in 1:nthreads
        for iλl in 1:nλl 
            int_f[iλl] += int_ft[tid][iλl]
        end
    end
end

@doc raw"""
    $I = I(0) + ϵ/κ (1 - \exp(-κ z))$
    $k << 1: I = I(0) + ϵ z$
    $k >> 1: I = I(0) + ϵ / κ$
"""
function integrate_intensity_over_Δs(Iλb::Vector{Float64}, κb::Vector{Float64}, ϵb::Vector{Float64},  Δs::Float64, par)
    κΔs_limit        = par.r.κΔs_limit
    omit_absorb_emit = par.r.omit_absorb_emit

    #plt.plot(Iλb)

    if omit_absorb_emit == :omit_none
        κbΔs = κb .* Δs
        exp_κΔs = exp.(-κbΔs)
        Iλb[:] = @. ifelse(κbΔs < κΔs_limit, 
                        Iλb .* exp_κΔs .+ ϵb.*Δs, 
                        Iλb .* exp_κΔs .+ ϵb ./ κb.*(1.0 .- exp_κΔs))
    elseif omit_absorb_emit == :omit_emission
        Iλb[:] = @. Iλb + ϵb*Δs
    else omit_absorb_emit == :emission_absorption
        κbΔs = κb .* Δs
        exp_κΔs = exp.(-κbΔs)
        Iλb[:] = @. Iλb * exp_κΔs
    end
end

function add_background()
#        #@time begin
#        # add background ?
#        if spec.par.background > 1.0e-10
#            iw = floor(Int64, ((ΔλL_mean + ΔλD_mean) / dλ * spec.par.Δλ_factor))
#            # compute movning average
#            ma_κ = moving_average5(spec.κ_c, iw*2)
#            ma_ϵ = moving_average5(spec.ϵ_c, iw*2)
#            # add background
#            for iλ in 1:nb_λ
#                spec.κ_c[iλ] += ma_κ[iλ] * spec.par.background
#                spec.ϵ_c[iλ] += ma_ϵ[iλ] * spec.par.background
#            end
#            # ning average once
#            if iN == 1 && iθ == 1
#                save_intensity_as_hdf5(joinpath(spec.par.out_dir, "moving_average_kappa"), spec.λ, ma_κ)
#            end
#        end
end

"""
species = :H2O
moleculardata = molec_data_dict[species]
"""
function get_species_line_data(par, spec, moleculardata; renew_hdf5=false)
    datfiles = get_data_files()
    if renew_hdf5
        hitran_to_hdf5(spec, datfiles[spec][:csv], datfiles[spec][:hdf5], datfiles[spec][:hdf5_compact], par.w.λmin, par.w.λmax, length(moleculardata.iso_a))
    end
    hdf5_path = datfiles[spec][:hdf5_compact]
    line_data = LineData(hdf5_path);
    line_data
end

function get_line_data(par, molec_data_dict; renew_hdf5=false)
    line_data_dict = Dict{Symbol,LineData}()
    for spec in par.m.species
        line_data_dict[spec] = get_species_line_data(par, spec, molec_data_dict[spec]; renew_hdf5=renew_hdf5)
    end
    line_data_dict
end
