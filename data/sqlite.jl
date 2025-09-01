using SQLite
using CSV
using DataFrames
using Printf

function qfiles_to_df(root)
    qfiles = filter(x -> x[1] == 'q', readdir(root))
    qfiles = [x[2] for x in sort([(parse(Int32, split(f[2:end], ".")[1]), f) for f in qfiles]) ]

    function sppa(x)
        sx = split(strip(x),  r"\s+")
        [parse(Float64, sx[1]), parse(Float64, sx[2])]
    end

    T = []
    Q = []
    for f in qfiles
        data = [sppa(x) for x in readlines(joinpath(root, f))][1:3500]
        M = reduce(hcat, data)'
        println(f, " ", size(M))
        push!(T, M[:,1])
        push!(Q, M[:,2])
    end
    A = Matrix{Float64}(undef, 3500, length(Q)+1)
    A[:,1] = T[1]
    for i in eachindex(Q)
        A[:,i+1] = Q[i]
    end
    colsyms = [Symbol(splitext(x)[1]) for x in qfiles]
    insert!(colsyms, 1, :T)
    println(length(colsyms))
    df0 = DataFrame(A, :auto)
    rename!(df0, colsyms)
    df0
end

function TQ_data_to_db!(db, QCO2, colsyms)
    df0 = DataFrame(QCO2, :auto)
    SQLite.load!(df0, db, "CO2_QT")
end

function molec_data_to_db!(db, molec)
    root = joinpath(@__DIR__, molec, molec*"_Q")

    df_qfiles   = qfiles_to_df(root)
    SQLite.load!(df_qfiles  , db, @sprintf("%s_TQ", molec)) 

    df_Isotopes = CSV.read(joinpath(root, @sprintf("%s_Isotopes.txt", molec)), DataFrame)
    SQLite.load!(df_Isotopes, db, @sprintf("%s_Iso", molec)) 

    if isfile(@sprintf("%s_q7-q122.txt",  molec))
        df_q7_q122  = CSV.read(joinpath(root, @sprintf("%s_q7-q122.txt",  molec)), DataFrame)
        SQLite.load!(df_q7_q122 , db, @sprintf("%s_qs", molec)) 
    end
end
molec = "H2O"
dbpath = "MolecularData.db"
db = SQLite.DB(dbpath)
molec_data_to_db!(db, molec)

CO2_rwfmt.csv
H2O_rwfmt.csv

df_lines = CSV.read(joinpath(root, @sprintf("%s_rwfmt.csv", molec)), DataFrame)
SQLite.load!(df_Isotopes, db, @sprintf("%s_linedata", molec)) 

df0_CO2_TQ = DataFrame(TQCO2, :auto)
rename!(df0_CO2_TQ, colsyms)
SQLite.load!(df0_CO2_TQ, db, "CO2_QT")df_qfiles  
df_Isotopes
df_q7_q122 


SQLite.load!(df_Isotopes, db, @sprintf("%Isotopes", molec))
SQLite.load!(df_q7_q122 , db, @sprintf("%q7_q122",  molec))



data[1]
reduce(vcat, data)
reduce(vcat,transpose.(x))
data[1]
split("5000       887489.52300000", r"\s+")
data[1]
db = SQLite.DB(dbpath)
SQLite.createtable!(db, "results", Tables.Schema(colnames, coltypes))
rowqm = "?"*repeat(", ?", length(columns.name)-1)


dbpath = "CO2.db"
