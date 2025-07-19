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

function ResultDB(dbpath, csih)
    if isfile(dbpath)
        @warne "database", dbpath, "exists"
    end

    colnames = ["hdf5_path", "ic", "iθ", "ih", "h", "θ",   "T",  "N", "species"]
    coltypes = [String, Int,  Int,   Int, Real, Real, Real,  Real, String]
    
    species = collect(keys(csih))
    for spec in species
        colnames = cat(colnames, ["cih$spec", "ΔλL$spec", "ΔλD$spec", "int_I$spec", "int_ϵ$spec", "int_κ$spec", "int_Iκ$spec"], dims=1)
        coltypes = cat(coltypes, [Real, Real, Real, Real, Real, Real, Real], dims=1)
    end
    
    #length(colnames)
    #length(coltypes)

    db = SQLite.DB(dbpath)
    SQLite.createtable!(db, "results", Tables.Schema(colnames, coltypes))

    columns = SQLite.columns(db, "results")
    #length(columns.name)

    io = IOBuffer()
    write(io, @sprintf("%s", columns.name[1]))
    for n in columns.name[2:end]
        write(io, @sprintf(", %s", n))
    end
    cln = String(take!(io))

    rowqm = "?"*repeat(", ?", length(columns.name)-1)

    ResultDB(dbpath, db, cln, rowqm, species)
end

function create_results_db(par)
    rdb = ResultDB(par[:paths][:dbpath], par[:c_ppm])
    rdb
end

function open_db(dbpath)
    SQLite.DB(dbpath)
end

function insert_into_rdb(rdb::ResultDB, ic, iθ, ih, h, θ, T, N, cihic, ΔλL_mean, ΔλD_mean, int_I, int_ϵ, int_κ, int_Iκ, htf5_path)
    species = collect(keys(cihic))
    row = Array{Any}([htf5_path, ic, iθ, ih, h, θ, T, N, list_to_string(rdb.species)])
    for spec in rdb.species
        row = cat(row, [cihic[spec], ΔλL_mean[spec], ΔλD_mean[spec], int_I[1], int_ϵ[1][spec], int_κ[1][spec], int_Iκ[1][spec]], dims=1)
    end

    DBInterface.execute(rdb.db, @sprintf("INSERT INTO results (%s) VALUES (%s);", rdb.colnames, rdb.rowqm), row)
end

function select_from_rdb(db; ic=1, iθ=1)
    result = DBInterface.execute(db, "SELECT * FROM results WHERE ic = $ic AND iθ = $iθ;")
    DataFrame(result)
end
