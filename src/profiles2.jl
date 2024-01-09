using RWPhysConst
using RWUtils
using RWFileIO
using StaticArrays

struct ExpPC
    n  :: Int64
    a  :: Float64
    x0 :: Float64
    x  :: Vector{Float64}
    y  :: Vector{Float64}
end
function ExpPC(n::Int64, x0)
    x = mgrid(-x0, x0, n)
    y = exp.(-x.^2)
    a = Float64(n-1) / (2.0*x0)
    ExpPC(n, a, -x0, x, y)
end
"""
exp(- c * (λ - λ0)^2) = exp(- ((λ - λ0)/d)^2)
c = 1/d^2
"""
function get_exp(e::ExpPC, λ, λ0, c)
    x = (λ-λ0) * c
    i = max(1, min(e.n, floor(Int64, (x - e.x0) * e.a) + 1))
    e[i]
end

struct GaussProfile
    λ0     :: Float64
    Δλ0    :: Float64
    sqrt_c :: Float64
    a      :: Float64
    e      :: ExpPC
end

function GaussProfile(M::Float64, T::Float64, λ0::Float64)
    ln2 = log(2.0)
    Δλ0 = sqrt(2.0 * c_kB * T / M) / c_c
    dπ  = sqrt(1.0/π)
    Δλλ = Δλ0 * λ0
    sqrt_c = sqrt(ln2) / Δλλ
    a = dπ * sqrt_c
    e = ExpPC(201, 5.0)
    GaussProfile(λ0, Δλ0, sqrt_c, a, e)
end

function gauss_profile(M::Float64, T::Float64, λ0::Float64)
    ln2 = log(2.0)
    Δλ0 = sqrt(2.0 * c_kB * T / M) / c_c
    dπ  = sqrt(1.0/π)
    Δλλ = Δλ0 * λ0
    sqrt_c = sqrt(ln2) / Δλλ
    a = dπ * sqrt_c
    SA_F64[λ0, Δλ0, sqrt_c, a]
end

@inline function fa(gp::GaussProfile, λ)
    x = ((λ - gp.λ0)*gp.sqrt_c)
    i = floor(Int64, (x-gp.e.x0) * gp.e.a + 1.0)
    return x,i
#    if i < 1
#        return 0.0
#    elseif i >= gp.e.n 
#        return 0.0
#    end
#    u = (x-gp.e.x[i]) * gp.e.a
#    gp.a * ((1.0-u)*gp.e.y[i] + u*(gp.e.y[i+1]))
end

@inline function fa(gp::GaussProfile, λ)
    x = ((λ - gp.λ0)*gp.sqrt_c)
    i = floor(Int64, (x-gp.e.x0) * gp.e.a + 1.0)
    if i < 1
        return 0.0
    elseif i >= gp.e.n 
        return 0.0
    end
    u = (x-gp.e.x[i]) * gp.e.a
    gp.a * ((1.0-u)*gp.e.y[i] + u*(gp.e.y[i+1]))
end

@inline function f(gp::GaussProfile, λ)
    @. gp.a * exp(- ((λ - gp.λ0)*gp.sqrt_c)^2)
end

function fv(g::GaussProfile, λ::Vector{Float64}, λ0)
    @. gp.a * exp(- ((λ - λ0)*gp.sqrt_c)^2)
end


struct LorentzProfile
    ΔλL :: Float64
    πΔλL :: Float64
end
function LorentzProfile(ΔλL)
    LorentzProfile(ΔλL, π * ΔλL)
end

function lorentz_profile(ΔλL)
    SA_F64[ΔλL, π*ΔλL]
end

function f(l::LorentzProfile, λ::Float64, λ0::Float64)
    dλ = (λ - λ0)/l.ΔλL
    1.0 / (l.πΔλL * (dλ^2 + 1.0))
end

@inline function fv(l::LorentzProfile, λ::Vector{Float64}, λ0)
    dλ = @. (λ - λ0)/l.ΔλL
    @. 1.0 / (π * l.ΔλL * (dλ^2 + 1.0))
end

mutable struct VoigtProfile
    Gauss   :: GaussProfile
    Lorentz :: LorentzProfile
end

"""
    Simple Voigt profile
"""
@inline function voigt(g::GaussProfile, l::LorentzProfile, λ, λ0)
    v = l.ΔλL / (g.Δλ0*g.λ0)
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return f(l, λ, λ0)
    end
    return v * f(l, λ, λ0) + (1.0 - v) * f(g, λ)
end

@inline function f(vp::VoigtProfile, λ, λ0)
    v = vp.Lorentz.ΔλL / (vp.Gauss.Δλ0*λ0)
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return f(vp.Lorentz, λ, λ0)
    else
        return v * f(vp.Lorentz, λ, λ0) + (1.0 - v) * f(vp.Gauss, λ, λ0)
    end
end

function fv(vp::VoigtProfile, λ, λ0)
    v = vp.Lorentz.ΔλL / (vp.Gauss.Δλ0*λ0)
    v = 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3;
    if v > 1.0
        return fv(vp.Lorentz, λ, λ0)
    else
        return v * f(vp.Lorentz, λ, λ0) + (1.0 - v) * f(self.Gauss, λ)
    end
end


function test_gauss()
    NO = 16.0
    NC = 12.0
    M = c_mp * NO * 2.0 + c_mp * NC
    T = c_TCK
    λ0 = 1.5e-5
    gp = GaussProfile(M, T, λ0)

    
    f1 = fa(gp, gp.λ0)
    f2 = f(gp, gp.λ0)
    println(f1, " ", f2)
end

test_gauss()

NO = 16.0
NC = 12.0
M = c_mp * NO * 2.0 + c_mp * NC
T = c_TCK
λ0 = 1.5e-5
gp = GaussProfile(M, T, λ0)
gp.a
λ = gp.λ0-gp.Δλ0*1.0e-5
λ = mgrid(gp.λ0-gp.Δλ0*1.0e-5, gp.λ0+gp.Δλ0*1.0e-5, 21)

@benchmark fa(gp, λ)
@benchmark f(gp,  λ)

for i in 1:21
    f1 = fa(gp, λ[i])
    f2 = f(gp,  λ[i])
    println(f1, " ", f2)
end

gp.e.y

lp = LorentzProfile(gp.Δλ0)

@benchmark voigt(gp, lp, λ, gp.Δλ0)

Δλ0 = gp.Δλ0
@benchmark f(lp, λ, Δλ0)


gp = gauss_profile(M, T, λ0)
function f(gp::SVector{4, Float64}, λ::Float64)
    #SA_F64[λ0, Δλ0, sqrt_c, a]
    gp[4] * exp(- ((λ - gp[1])*gp[3])^2)
end
@benchmark f(gp, λ)

function fL(lp, λ, λ0)
    dλ = (λ - λ0)/lp.ΔλL
    1.0 / (lp.πΔλL * (dλ^2 + 1.0))
end
@benchmark fL(lp, λ, λ0)
