using PhysConst
using SimpleLog
using StaticArrays

function GaussProfile(M::Float64, T::Float64, λ0::Float64)
    ln2 = log(2.0)
    Δλ0 = sqrt(2.0 * c_kB * T / M) / c_c
    dπ  = sqrt(1.0/π)
    Δλλ = Δλ0 * λ0
    sqrt_c = sqrt(ln2) / Δλλ
    a = dπ * sqrt_c
    SA_F64[λ0, Δλ0, sqrt_c, a]
end

@inline function fG(gp, λ)
    gp[4] * exp(- ((λ - gp[1])*gp[3])^2)
end

function LorentzProfile(ΔλL)
    SA_F64[ΔλL, π*ΔλL]
end

function fL(lp, λ, λ0)
    dλ = (λ - λ0)/lp[1]
    1.0 / (lp[2] * (dλ^2 + 1.0))
end

"""
    Simple Voigt profile
"""
@inline function voigt(gp, lp, λ, λ0)
    v = lp[1] / (gp[2] * gp[1])
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return fL(lp, λ, λ0)
    end
    return v * fL(lp, λ, λ0) + (1.0 - v) * fG(gp, λ)
end
