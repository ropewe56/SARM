
struct GaussProfile
    Δλ0  :: Float64
    ln2  :: Float64
    dπ   :: Float64
end

function GaussProfile(M::Float64, T::Float64)
    ln2 = log(2.0)
    Δλ0 = sqrt(2.0 * c_kB * T / M) / c_c
    dπ  = sqrt(1.0/π)
    GaussProfile(Δλ0, ln2, dπ)
end

@inline function f_gauss(gp::GaussProfile, λ, λ0)
    @fastmath Δλ = gp.Δλ0 * λ0
    @fastmath c  = gp.ln2 / Δλ^2
    @fastmath gp.dπ * sqrt(c) * exp(- c * (λ - λ0)^2)

end

struct LorentzProfile
    ΔλL :: Float64
end

@inline function f_lorentz(l::LorentzProfile, λ, λ0)
    @fastmath dλ = (λ - λ0)/l.ΔλL
    @fastmath 1.0 / (π * l.ΔλL * (dλ^2 + 1.0))
end


mutable struct VoigtProfile
    Gauss   :: GaussProfile
    Lorentz :: LorentzProfile
end

"""
    Simple Voigt profile
"""
@inline function voigt(ΔλG, ΔλL, fg, fl, λ, λ0)
    v = ΔλL / ΔλG
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return fl
    end
    @. v * fl + (1.0 - v) * fg
end

@inline function fV(vp::VoigtProfile, λ, λ0)
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
        return v * fa(vp.Lorentz, λ, λ0) + (1.0 - v) * fa(self.Gauss, λ, λ0)
    end
end

