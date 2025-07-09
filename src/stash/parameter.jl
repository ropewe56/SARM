using JSON
using Parameters

struct MoleculeData
    isotope_names       :: Vector{String}
    isotope_masses      :: Vector{Float64}
    TQ                  :: Matrix{Float64}
    ids                 :: Matrix{Int64}
    spectral_data       :: Matrix{Float64}
    max_isotope_ids     :: Vector{Int64}
    CO2_c0_ppm          :: Float64 
    H2O_c0_ppm          :: Float64 
end

@with_kw mutable struct RunParameter
    λmin                :: Float64         = 14.0e-6
    λmax                :: Float64         = 16.0e-6
    Δλ                  :: Float64         = 1.0e-11
    ΔλL                 :: Float64         = 1.0e-11
    Δλ_factor           :: Float64         = 1.0
    surface_T           :: Float64         = 300.0
    Planck_Ts           :: Vector{Float64} = zeros(Float64, 0)
    initial_intensity   :: String          = "planck"
    p_ref               :: Float64         = 1.06e5
    T_ref               :: Float64         = 300.0
    albedo              :: Float64         = 0.3
    background          :: Float64         = 1.0
    T_of_h              :: Bool            = true
    N_of_h              :: Bool            = true
    with_emission       :: Bool            = true
    integrate           :: Bool            = true
    out_dir             :: String          = ""
    logfile             :: String          = ""
    infofile            :: String          = ""
    z_iout              :: String          = ""
    CO2it_ppm           :: Matrix{Float64}([1.0, 2.0, 3.0])
    θ_deg               :: Vector{Float64}([0.0, 40.0, 80.0])
    z                   :: Dict("zmin"=> 0.0, "zmax"=> 7.0e4, "dzmin"=> 10.0, "dzmax"=> 1.0e3, "n" => 200,  "exponent" => 1)
end

function Parameter(parameter_json_file)
    # 
    json = JSON.parsefile(parameter_json_file; dicttype=Dict, inttype=Int64, use_mmap=true)
    mCO2                = json["run_parameter"]["CO2_mass" ]
    λmin                = json["run_parameter"]["λmin"]
    λmax                = json["run_parameter"]["λmax"]
    Δλ                  = json["run_parameter"]["Δλ"]
    ΔλL                 = json["run_parameter"]["ΔλL"]
    Δλ_factor           = json["run_parameter"]["Δλ_factor"]
    surface_T           = json["run_parameter"]["surface_T"]
    Planck_Ts           = json["run_parameter"]["Planck_Ts"]
    initial_intensity   = json["run_parameter"]["initial_intensity"]
    T_ref               = json["run_parameter"]["T_ref"]
    p_ref               = json["run_parameter"]["p_ref"]
    albedo              = json["run_parameter"]["albedo"]
    T_of_h              = json["run_parameter"]["T_of_h"]
    N_of_h              = json["run_parameter"]["N_of_h"]
    with_emission       = json["run_parameter"]["with_emission"]
    background          = json["run_parameter"]["back_ground"]
    max_isotope_id      = json["run_parameter"]["max_isotope_id"]
    integrate           = json["run_parameter"]["integrate"]
    out_dir             = json["paths"]["out_dir"]
    logfile             = json["paths"]["logfile"]
    infofile            = json["paths"]["infofile"]
    NCO2_path           = json["paths"]["NCO2_path"]
    theta_path          = json["paths"]["theta_path"]
    specdat_path        = json["paths"]["specdat_path"]
    T_Q_path            = json["paths"]["T_Q_path"]
    h_p_path            = json["paths"]["h_p_path"]
    h_T_path            = json["paths"]["h_T_path"]
    z_path              = json["paths"]["z_path"]
    z_iout              = json["paths"]["z_iout"]

    mkpath(joinpath(out_dir, "intensity"))
    mkpath(joinpath(out_dir, "spectrum"))

    Parameter(  mCO2              ,
                λmin              ,
                λmax              ,
                Δλ                ,
                ΔλL               ,
                surface_T         ,
                Planck_Ts         ,
                initial_intensity ,
                T_of_h            ,
                N_of_h            ,
                Δλ_factor         ,
                p_ref             ,
                T_ref             ,
                albedo            ,
                with_emission     ,
                background        ,
                max_isotope_id    ,
                integrate         ,
                out_dir           ,
                logfile           ,
                infofile          ,
                NCO2_path         ,
                theta_path        ,
                specdat_path      ,
                T_Q_path          ,
                h_p_path          ,
                h_T_path          ,
                z_path            ,
                z_iout            )
end

mutable struct Hitran
    width :: Vector{Int64}
    pos   :: Vector{Int64}
end

function Hitran(par::Parameter)
    width1 = [ 2,  1,  12,  10,  10,   5,   5,  10,   4,   8,  15,  15,  15,  15]
    width2 = [ 1,  2,   1,   7,   7]
    w1 = size(width1, 1)
    w2 = size(width2, 1)
    n = 160 - w1 - w2
    width = [  2,   1,  12,  10,  10,   5,   5,  10,   4,   8,  15,  15,  15,  15, n,  1,  2,  2,  7,  7]

    pos = Vector{Int64}(undef,0)
    sum_of_elems = 0.0
    for n in 1:eachindex(width)
        push!(pos, sum_of_elems)
        sum_of_elems += width[n]
    end
    pos[end] = pos[end] - 1
    Hitran(width, pos)
end

# molec_id,local_iso_id,nu,sw,a,gamma_air,gamma_self,elower,n_air,delta_air,gp,gpp

function read_hitran_data(hit::Hitran, filename)
    lines = reaΔλines(filename)

    hitran_data = Vector{Vector{Float64}}(undef,0)

    for line in lines[2:end]
        data = split(line, ",")
        if size(data,1) == 12
            mid    = parse(Float64, data[1])
            iid    = parse(Float64, data[2])
            ν      = parse(Float64, data[3])
            S      = parse(Float64, data[4])
            A      = parse(Float64, data[5])
            γ_a    = parse(Float64, data[6])
            γ_s    = parse(Float64, data[7])
            ν_l    = parse(Float64, data[8])
            n_a    = parse(Float64, data[9])
            δ_a    = parse(Float64, data[10])
            gp     = parse(Float64, data[11])
            gpp    = parse(Float64, data[12])

            ν      = ν   / 1.0e-2 # 1/m
            ν_l    = ν_l / 1.0e-2 # 1/m
            γ_a    = γ_a / 1.0e-2 # 1/m
            γ_s    = γ_s / 1.0e-2 # 1/m

            λ      = 1.0 / ν
            ΔE_ul  = hit.h * hit.c / λ
            E_l    = hit.h * hit.c * ν_l
            E_u    = E_l + ΔE_ul
            Δλa    = γ_a * λ^2 # hwh
            Δλs    = γ_s * λ^2 # hwh

            # S = 1/N * (N_l * B_lu - N_u * B_ul) * h * ν / c
            # A_ul = 8 * π * h * ν**3 * B_ul
            # g_l * B_lu = g_u * B_ul
            # S_ν^N(T) = 1/Q_{tot}(T) * ( exp(-E_l/(kB*T)) * gl * B_lu - exp(-E_u/(kB*T)) * gu * B_ul) * h * ν_0 / c
            # S_ν^N(T) = g_1/Q_{tot}(T) * A_ul / (8 π c \nu_0**2) exp(-E_l/(kB*T)) * (1.0 - epx(-h ν / (kB*T)))

            c = zeros(Float64, 9)
            c[ 1] = λ                      # m
            c[ 2] = E_l                    # eV
            c[ 3] = E_u                    # eV
            c[ 4] = S                      # 1/s
            c[ 5] = A                      # 1/s
            c[ 6] = γ_a * 1.0e-5           # 1/m/pascal
            c[ 7] = γ_s * 1.0e-5           # 1/m/pascal
            c[ 8] = n_a                    #
            c[ 9] = δ_a / 1.0e-2 * 1.0e-5  #  δ_a : [1/cm 1/atm] => [1/1e-2m 1/10e5p] => 1/m * 1/pascal
            push!(hitran_data, c)
        end
    end
    hitran_data
end
