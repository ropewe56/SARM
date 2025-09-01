using SQLite
using CSV
using DataFrames
using Printf
using PhysConst

function to_MKS(df)
    hc    = c_h * c_c
    cm    = 1.0e-2
    atm   = 1.01325e5

    iso_id = df[!,:local_iso_id]
    λ21   = cm ./ (df[!,"nu"])
    S21r  = df[!,:sw] .* cm 
    A     = df[!,:a]
    γair  = df[!,:gamma_air]  ./ (cm * atm)
    γself = df[!,:gamma_self] ./ (cm * atm)
    E1    = (hc * cm) .* (df[!,:elower])
    nair  = df[!,:n_air]
    δair  = df[!,:delta_air]  ./ (cm * atm)
    g2    = df[!,:gp]
    g1    = df[!,:gpp]

    colnames = ["iso_id", "λ21", "S21r", "A21", "γair", "γself", "E1", "nair", "δair", "g2", "g1"]
    values =   [ iso_id, 
                 λ21, 
                 S21r,
                 A,
                 γair ,
                 γself,
                 E1   ,
                 nair ,
                 δair ,
                 g2   ,
                 g1   ]

    df2 = DataFrame(colnames .=> values)
    sort!(df2, [:λ21]);
    df2
end

function create_unique_index(db, table_name, col_names)
    cmd = @sprintf("CREATE UNIQUE INDEX unique_idx ON %s (%s);", table_name, join(" ,", col_names))
    DBInterface.execute(db, cmd)
end

function qfiles_to_df(molec)
    root = joinpath(@__DIR__, molec, molec*"_Q")
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
    TQ = Matrix{Float64}(undef, 3500, length(Q)+1)
    TQ[:,1] = T[1]
    for i in eachindex(Q)
        TQ[:,i+1] = Q[i]
    end

    colsyms = [Symbol(@sprintf("Q_%0d", i)) for i in eachindex(qfiles)]
    insert!(colsyms, 1, :T)
    df0 = DataFrame(TQ, colsyms)
    df0
end

function molec_data_to_db!(db, molec)
    root = joinpath(@__DIR__, molec, molec*"_Q")
    mass_factor = 1.0e-3/6.02214076e23

    df_Q = qfiles_to_df(molec)

    df_iso = if isfile(joinpath(root, @sprintf("%s_q7-q122.txt",  molec)))
        df_qq = CSV.read(joinpath(root, @sprintf("%s_q7-q122.txt",  molec)), DataFrame)
        miso = df_qq[:,Symbol("MolarMass/g·mol-1")] .* mass_factor
        aiso = df_qq.Abundance 
        iso_id = df_qq.localID
        DataFrame(iso_id = iso_id, aiso = aiso, misos = miso)
    else
        df = CSV.read(joinpath(root, @sprintf("%s_Isotopes.txt", molec)), DataFrame)
        DataFrame(iso_id = collect(1:nrow(df)), aiso = df.IsoAbundance, miso=df[!,Symbol("MolarMass(g)")] .* mass_factor)
    end

    table_name = @sprintf("%s_TQ", molec)
    SQLite.load!(df_Q, db, table_name)
    cmd = @sprintf("CREATE UNIQUE INDEX unique_%s ON %s (%s);", table_name, table_name, "T")
    DBInterface.execute(db, cmd)

    table_name = @sprintf("%s_Iso", molec)
    SQLite.load!(df_iso, db, table_name) 
    cmd = @sprintf("CREATE UNIQUE INDEX unique_%s ON %s (%s);", table_name, table_name, "iso_id")
    DBInterface.execute(db, cmd)
end

function is_unique(df)
    nr = nrow(df)
    nc = ncol(df)    
    for i in 1:nr-1
        flag = true
        for j in 1:nc
            if df[i,j] != df[i+1,j]
                flag = false
            end
        end
        if flag println(i, "\n ",df[i,:], df[i+1,:]) end
    end
end

function line_data_db!(db, molec)
    root = joinpath(@__DIR__, molec)
    df = CSV.read(joinpath(root, @sprintf("%s_rwfmt.csv", molec)), DataFrame)
    unique!(df)
    
    df2 = to_MKS(df)

    table_name = @sprintf("%s_linedata", molec)
    schema = Tables.Schema(names(df2), eltype.(eachcol(df2)))    
    SQLite.createtable!(db, table_name, schema; temp=false, ifnotexists=true)
    cmd = @sprintf("CREATE UNIQUE INDEX unique_%s ON %s (%s);", table_name, table_name, join(names(df2), ", "))
    DBInterface.execute(db, cmd)

    SQLite.load!(df2, db, @sprintf("%s_linedata", molec)) 
end

function raw_data_to_db()
    dbpath = joinpath(@__DIR__, "MolecularData.db")
    db = SQLite.DB(dbpath)
    for molec in ["CO2", "H2O"]
        molec_data_to_db!(db, molec)
        line_data_db!(db, molec)
    end
end
#raw_data_to_db()
