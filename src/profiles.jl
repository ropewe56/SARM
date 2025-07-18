using PhysConst
using SimpleLog

const LOG2 = log(2.0)

@inline function adapt_f(λ, λ0, Δλh, f, f_adapt)
    if f_adapt == :none
        return f
    end
    
    df = min(f[1], f[end])
    if f_adapt == :scaletail
        (f .- df) .* ( 1.0 / (1.0 - df*length(f)/sum(f)) * 0.5*π/max(atan((λ0 - λ[1])/Δλh), atan((λ[end] - λ0)/Δλh)) )
    elseif  f_adapt == :tail
        # tail_energy
        f .* 0.5*π/max(atan((λ0 - λ[1])/Δλh), atan((λ[end] - λ0)/Δλh))
    else
        # scale
        (f .- df) ./ (1.0 - df*length(f)/sum(f))
    end
end

@inline function f_gauss(λ::Vector{Float64}, λ0, ΔλGh, fG_adapt)
    a = LOG2/ΔλGh^2
    f = @. sqrt(a/π) * exp(- a * (λ - λ0)^2)
    adapt_f(λ, λ0, ΔλGh, f, fG_adapt)
end

@inline function f_lorentz(λ::Vector{Float64}, λ0, ΔλLh, fL_adapt)
    f = @. 1.0 / (π * ΔλLh * (1.0 + ((λ - λ0)/ΔλLh)^2))
    adapt_f(λ, λ0, ΔλLh, f, fL_adapt)
end

@inline function voigt(λ::Vector{Float64}, λ0, ΔλLh, ΔλGh, fL_adapt, fG_adapt)
    fL = f_lorentz(λ, λ0, ΔλLh, fL_adapt)
    v = ΔλLh / ΔλGh
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return fL
    end
    fG = f_gauss(λ, λ0, ΔλGh, fG_adapt)
    @. v * fL + (1.0 - v) * fG
end

