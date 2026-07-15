using CSV
using DataFrames
using PhysConst
using PhysConst.UnitConst
using SQLite
using Statistics

function convert_units(row)
    _cm   = 1.0e-2
    _atm  = 1.01325e5
    hc    = c_h*c_c
    NCO20 = 430.0e-6*2.5e25

    molid = row[:molec_id]  
    isoid = row[:local_iso_id]

    nu21  = row[:nu]
    A21   = row[:a]
    γair  = row[:gamma_air]
    γself = row[:gamma_self]
    elower= row[:elower]
    nair  = row[:n_air]
    δair  = row[:delta_air]
    g2    = row[:gp]
    g1    = row[:gpp]

    nu21_m = nu21 ./ _cm                                               # [1/m]
    λ210   = 1.0 ./ nu21_m                                             # [m]
    E1     = hc * elower ./ _cm                                               # [1/m]
    ΔE21   = hc * nu21_m
    E2     = E1 + ΔE21

    B21 = A21 * λ210^3 / (8π * c_h)
    B12 = g2 / g1 * B21
    σ   = c_h * λ210/c_c * B12 / λ210  # Js/m m^3/Js/s^2
    κ   = c_h * λ210/c_c * B12
    ϵ   = c_h * c_c/λ210 * A21
    Λ   = 1/(NCO20*σ)

    γair  = γair  ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal
    γself = γself ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal
    δair  = δair  ./ (_cm * _atm)                                    # [1/(m*Pa)] , cm => m, atm => pascal

    [molid, isoid, λ210, ΔE21, E1, E2, A21, B21, B12, g2, g1, κ, ϵ, σ, Λ,
        γair, γself, nair, δair]
end

function transfrom_data(in_path)
    CO2 = CSV.read(in_path, DataFrame)

    df = DataFrame(molid=Int32[], isoid=Int32[], λ210=Float64[], ΔE21=Float64[], E1=Float64[], E2=Float64[], 
        A21=Float64[], B21=Float64[], B12=Float64[], g2=Float64[], g1=Float64[], 
        κ=Float64[], ϵ=Float64[], σ=Float64[], Λ=Float64[],
        γair=Float64[], γself=Float64[], nair=Float64[], δair=Float64[])

    for row in eachrow(CO2)
        push!(df, convert_units(row))
    end

    sort!(df, [:λ210])
    df
end

function load_all_line_data()
    dbpath = joinpath(@__DIR__, "line_data.sqlite3")
    rm(dbpath)

    db = SQLite.DB(dbpath)

    mol = ("CO2", "CH4", "H2O", "N2O")
    for m in mol
        in_path = joinpath(@__DIR__, m*".csv")
        table_name = m*"_linedata"
        df = transfrom_data(in_path)

        SQLite.load!(df, db, table_name)
    end
end

function get_mfp()
    dbpath = joinpath(@__DIR__, "line_data.sqlite3")
    db = SQLite.DB(dbpath)
    table_name = "CO2_linedata"

    df = DBInterface.execute(db, "SELECT * FROM $table_name") |> DataFrame

    Λ = df[!,:Λ]
    @info extrema(Λ), median(Λ)

    E1 = df[!,:E1]
    @info extrema(E1), median(E1)

    E2 = df[!,:E2]
    @info extrema(E2), median(E2)

    g1 = df[!,:g1]
    N0 = 430.0e-6*2.5e25
    T = 300.0
    c_kB*T/c_e
    N1 = @. exp(-E1 / (c_kB*T)) * N0 * g1
    @info extrema(N1), median(N1)
    sum(@. ifelse(N1 > 1.0e18, 1, 0))
end

function check_units(db, molec)
    table_name = molec*"_linedata"
    df = DBInterface.execute(db, "SELECT * FROM $table_name") |> DataFrame
    names(df)
    row = df[1,:]

    λ = 13.0e-4
    λ = 17.0e-4
    nu_min = 780.0
    nu_max = 580.0

    N0 = 430.0e-6*2.5e25 / u_m^3

    molid = row[:molid]  
    isoid = row[:isoid]

    A21 = row[:A21]/u_s
    λ   = row[:λ210]  * u_m
    B21 = A21 * λ^3/(8π*u_h)
    B12 = row[:g2]/row[:g1] * B21
    σ   = u_h * λ/u_c * B12 / λ  # Js/m m^3/Js/s^2
    κ   = u_h * λ/u_c * B12
    ϵ   = u_h * u_c/λ * A21
    Λ   = 1/(N0*σ)

    γair  = row[:γair]
    γself = row[:γself]
    E1    = row[:E1]
    nair  = row[:nair]
    δair  = row[:δair]
    g2    = row[:g2]
    g1    = row[:g1]

end