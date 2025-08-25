using StaticArrays

mutable struct Profile
    nf     :: Int64
    nth    :: Int64
    ft     :: Vector{Vector{Float64}}
end

function make_fs(nf, nthreads, nλl)
    ft = Vector{Vector{Float64}}(undef, nthreads)
    for i in 1:nthreads
        ft[i] = zeros(Float64, nf)
    end
    ft
end

function Profile(nf, nthreads, nλl)
    ft = make_fs(nf, nthreads, nλl)
    Profile(nf, nthreads, ft)
end

function update_profile(p::Profile, nf, nthreads, nλl)
    if nf > p.nf || nthreads > p.nth
        p.ft = make_fs(nf, nthreads, nλl)
    end
end

mutable struct PreAlloc
    linedata_dict :: Dict{Symbol,Vector{SVector{13,Float64}}}
    ϵb            :: Vector{Float64}
    κb            :: Vector{Float64}
    ϵbt           :: Vector{Vector{Float64}}
    κbt           :: Vector{Vector{Float64}}
    ϵbs           :: Dict{Symbol, Vector{Float64}}
    κbs           :: Dict{Symbol, Vector{Float64}}
    int_fs        :: Dict{Symbol, Vector{Float64}}
    int_fst       :: Dict{Symbol, Vector{Vector{Float64}}}
    pr            :: Profile
end

function PreAlloc(par, line_data_dict, nλb)

    nthreads = Threads.nthreads()

    ϵb  = zeros(Float64, nλb)
    κb  = zeros(Float64, nλb)

    ϵbt = Vector{Vector{Float64}}(undef,nthreads)
    κbt = Vector{Vector{Float64}}(undef,nthreads)
    for i in 1:nthreads
       ϵbt[i] = zeros(Float64, nλb)
       κbt[i] = zeros(Float64, nλb)
    end

    linedata_dict = Dict{Symbol,Vector{SVector{13,Float64}}}()
    int_fs = Dict{Symbol, Vector{Float64}}()
    int_fst = Dict{Symbol, Vector{Vector{Float64}}}()
    ϵbs = Dict{Symbol, Vector{Float64}}()
    κbs = Dict{Symbol, Vector{Float64}}()

    nλl_max = 0
    for spec in par.m.species
        nλl                 = get_nλl(line_data_dict[spec])
        nλl_max             = max(nλl_max, nλl)
        linedata_dict[spec] = Vector{SVector{13,Float64}}(undef, nλl)
        int_fs[spec]        = zeros(Float64, nλl)
        int_fst[spec]       = Vector{Vector{Float64}}(undef, nthreads)
        for i in 1:nthreads
            int_fst[spec][i] = zeros(Float64, nλl)
        end
        ϵbs[spec]           = zeros(Float64, nλb)
        κbs[spec]           = zeros(Float64, nλb)
    end

    pr = Profile(1000, nthreads, nλl_max)

    PreAlloc(linedata_dict,
             ϵb           ,
             κb           ,
             ϵbt          ,
             κbt          ,
             ϵbs          ,
             κbs          ,
             int_fs       ,
             int_fst      ,
             pr)
end

