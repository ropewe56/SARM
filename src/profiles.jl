using PhysConst
using SimpleLog
using LoopVectorization

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

@inline function f_gauss(λ, λ0, ΔλGh, fG_adapt)
    a = LOG2/ΔλGh^2
    f = Vector{Float64}(undef, length(λ))
    @turbo for i in eachindex(λ)
        f[i] = sqrt(a/π) * exp(- a * (λ[i] - λ0)^2)
    end
    #adapt_f(λ, λ0, ΔλGh, f, fG_adapt)
    f
end

@inline function f_lorentz(λ, λ0, ΔλLh, fL_adapt)
    f = Vector{Float64}(undef, length(λ))
    @turbo for i in eachindex(λ)
        f[i] = 1.0 / (π * ΔλLh * (1.0 + ((λ[i] - λ0)/ΔλLh)^2))
    end
    #adapt_f(λ, λ0, ΔλLh, f, fL_adapt)
    f
end

@inline function voigt(λ, λ0, ΔλLh, ΔλGh, fL_adapt, fG_adapt)
    fL = f_lorentz(λ, λ0, ΔλLh, fL_adapt)
    v = ΔλLh / ΔλGh
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
    if v > 1.0
        return fL
    end
    fG = f_gauss(λ, λ0, ΔλGh, fG_adapt)
    @. v * fL + (1.0 - v) * fG
end

