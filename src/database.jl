using SQLite
using DataFrames

mutable struct ResultDB
    dbpath
    db
    colnames
    rowqm
    species
end

function list_to_string(lst)
    io = IOBuffer()
    write(io, @sprintf("%s", lst[1]))
    for n in lst[2:end]
        write(io, @sprintf(", %s", n))
    end
    String(take!(io))
end

function ResultDB(par)
    dbpath = par[:paths][:dbpath]

    if isfile(dbpath)
        @warne @sprintf("Database %s exists. Will be deleted!", dbpath)
        rm(dbpath; force=true)
        #run(`rm $dbpath`)
    end
    db = SQLite.DB(dbpath)

    colnames = ["hdf5_path", "ic", "iθ", "ih", "h", "θ",  "T",  "N", "cCO2", "int_I", "int_ϵ", "int_Iκ", "species"]
    coltypes = [String, Int,  Int,   Int, Real, Real, Real,  Real, Real, Real, Real, Real,  Real, String]
    
    species = par[:species]
    for spec in species
        colnames = cat(colnames, ["cih$spec", "ΔλL$spec", "ΔλD$spec", "int_ϵs$spec", "int_Iκs$spec", "mean_κs$spec"], dims=1)
        coltypes = cat(coltypes, [Real, Real, Real, Real, Real, Real], dims=1)
    end

    SQLite.createtable!(db, "results", Tables.Schema(colnames, coltypes))

    columns = SQLite.columns(db, "results")

    @assert(length(colnames) == length(columns.name))

    io = IOBuffer()
    write(io, @sprintf("%s", columns.name[1]))
    for n in columns.name[2:end]
        write(io, @sprintf(", %s", n))
    end
    cln = String(take!(io))

    rowqm = "?"*repeat(", ?", length(columns.name)-1)

    ResultDB(dbpath, db, cln, rowqm, species)
end

function open_db(dbpath)
    SQLite.DB(dbpath)
end

function insert_into_resultdb(result_db::ResultDB, hdf5_path::String, ic::Int64, iθ::Int64, ih::Int64, 
                                h::Float64, θ::Float64, T::Float64, N::Float64, cihic, ΔλL_mean, ΔλD_mean, 
                                int_I, int_ϵ, int_Iκ, int_ϵs, mean_κs, int_Iκs)
    species = collect(keys(cihic))
    cCO2 = if haskey(cihic, :CO2)
        cihic[:CO2]
    else
        cihic[species[1]]
    end

    row = Array{Any}([hdf5_path, ic, iθ, ih, h, θ, T, N, cCO2, int_I[1], int_ϵ[1], int_Iκ[1], list_to_string(result_db.species)])
    for spec in result_db.species
        #@infoe [cihic[spec], ΔλL_mean[spec], ΔλD_mean[spec], int_ϵs[spec][1], int_Iκs[spec][1], mean_κs[spec][1]]
        row = cat(row, [cihic[spec], ΔλL_mean[spec], ΔλD_mean[spec], int_ϵs[spec][1], int_Iκs[spec][1], mean_κs[spec][1]], dims=1)
    end

    DBInterface.execute(result_db.db, @sprintf("INSERT INTO results (%s) VALUES (%s);", result_db.colnames, result_db.rowqm), row)
end

function select_from_rdb(db; ic=1, iθ=1)
    result = DBInterface.execute(db, "SELECT * FROM results WHERE ic = $ic AND iθ = $iθ;")
    DataFrame(result)
end
