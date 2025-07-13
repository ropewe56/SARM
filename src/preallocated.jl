using StaticArrays

mutable struct Preallocated
    κbt          :: Matrix{Float64}
    ϵbt          :: Matrix{Float64}
    fbt          :: Matrix{Float64}
    κb           :: Vector{Float64}
    ϵb           :: Vector{Float64}
    ΔλLs         :: Vector{Float64}
    ΔλDs         :: Vector{Float64}
    ΔλL_mean     :: Vector{Float64}
    ΔλD_mean     :: Vector{Float64}
    ld_pTN       :: Vector{SVector{9, Float64}}
    linedata_pTN :: Vector{Vector{SVector{9, Float64}}}
end
function Preallocated()
    κbt          = Matrix{Float64}(undef,0,0)
    ϵbt          = Matrix{Float64}(undef,0,0)
    fbt          = Matrix{Float64}(undef,0,0)
    κb           = Vector{Float64}(undef,0)
    ϵb           = Vector{Float64}(undef,0)
    ΔλLs         = Vector{Float64}(undef,0)
    ΔλDs         = Vector{Float64}(undef,0)
    ΔλL_mean     = Vector{Float64}(undef,0)
    ΔλD_mean     = Vector{Float64}(undef,0)
    ld_pTN       = Vector{SVector{9, Float64}}(undef,0)
    linedata_pTN = Vector{Vector{SVector{9, Float64}}}(undef,0)
    
    Preallocated(   κbt         ,
                    ϵbt         ,
                    fbt         ,
                    κb          ,
                    ϵb          ,
                    ΔλLs        ,
                    ΔλDs        ,
                    ΔλL_mean    ,
                    ΔλD_mean    ,
                    ld_pTN      ,
                    linedata_pTN)
end

@inline function a2(A, n1, n2, flag)
    A = if size(A) == (n1,n2) 
        A
    else
        Matrix{Float64}(undef, n1,n2) 
    end
    if flag
        fill!(A, 0.0)
    end
    A
end

@inline function alloc2(pre, name, n1, n2, flag=false)
    A = if name == :κbt          
        a2(pre.κbt, n1, n2, flag)
    elseif name == :ϵbt          
        a2(pre.ϵbt, n1, n2, flag)
    elseif name == :fbt          
        a2(pre.fbt, n1, n2, flag)
    end
    A
end

@inline function a1(A, n1, flag)
    A = if size(A,1) == n1
        A
    else
        Vector{Float64}(undef, n1) 
    end
    if flag
        fill!(A, 0.0)
    end
    A
end

@inline function alloc1(pre, name, n1, flag=false)
    A = if name == :ΔλLs    
        a1(pre.ΔλLs, n1, flag)
    elseif name == :ΔλDs
        a1(pre.ΔλDs, n1, flag)
    elseif name == :ΔλL_mean
        a1(pre.ΔλL_mean, n1, flag)
    elseif name == :ΔλD_mean
        a1(pre.ΔλD_mean, n1, flag)
    elseif name == :κb
        a1(pre.κb, n1, flag)
    elseif name == :ϵb
        a1(pre.ϵb, n1, flag)
    end
    A 
end

function alloc_ld_pTN(n)
    if length(ld_pTN) == n
        ld_pTN
    else
        Vector{SVector{9, Float64}}(undef,0)
    end
end

function linedata_pTN(n)
    if length(linedata_pTN) == n
        linedata_pTN
    else
        Vector{Vector{SVector{9, Float64}}}(undef,0)
    end
end
