using PhysConst
using SimpleLog
using LoopVectorization
using Bumper

@inline function adapt_f!(f, λ)
    imax = argmax(f)
    Δλ = λ[2] - λ[1]
    if f[1] < f[end]
        @. f[:] = f .- f[1]
        fsum = sum(f[1:imax])*Δλ
        @. f[:] = f[:] * 0.5/fsum
    else
        @. f[:] = f .- f[end]
        fsum = sum(f[imax:end])*Δλ
        @. f[:] = f[:] * 0.5/fsum
    end
end

#@inline function f_gauss(λ, λ0, ΔλGh)
#    a = LOG2/ΔλGh^2
#    f = @alloc(Float64, length(λ)) 
#    @turbo for i in eachindex(λ)
#        f[i] = sqrt(a/π) * exp(- a * (λ[i] - λ0)^2)
#    end
#    f
#end
#
#@inline function f_lorentz(λ, λ0, ΔλLh)
#    f = @alloc(Float64, length(λ)) 
#    @turbo for i in eachindex(λ)
#        f[i] = 1.0 / (π * ΔλLh * (1.0 + ((λ[i] - λ0)/ΔλLh)^2))
#    end
#    f
#end
#
#@inline function voigt(λ, λ0, ΔλLh, ΔλGh, f_adapt)
#    fL = f_lorentz(λ, λ0, ΔλLh)
#    v = ΔλLh / ΔλGh
#    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)
#    f = if v > 1.0
#        fL
#    else
#        fG = f_gauss(λ, λ0, ΔλGh)
#        @. v * fL + (1.0 - v) * fG
#    end
#    if f_adapt
#        adapt_f!(f, λ)
#    end
#end

@inline function f_lorentz!(f, λ, λ0, ΔλLh)
    @turbo for i in eachindex(λ)
        f[i] = 1.0 / (π * ΔλLh * (1.0 + ((λ[i] - λ0)/ΔλLh)^2))
    end
end

@inline function f_gauss!(f, v, λ, λ0, ΔλGh)
    a = LOG2/ΔλGh^2
    vv = (1.0 -  v)
    @turbo for i in eachindex(λ)
        f[i] = f[i] * v + sqrt(a/π) * exp(- a * (λ[i] - λ0)^2) * vv
    end
end

@inline function voigt!(f, λ, λ0, ΔλLh, ΔλGh, f_adapt)
    v = ΔλLh / ΔλGh
    v = max(0.0, 1.36606 * v - 0.47719 *v^2 + 0.11116 * v^3)

    f_lorentz!(f, λ, λ0, ΔλLh)
    if v <= 1.0
        f_gauss!(f, v, λ, λ0, ΔλGh)
    end

    #if f_adapt
    #    adapt_f!(f, λ)
    #end
end
