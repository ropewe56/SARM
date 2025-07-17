    # S f(ν) = σ
    # f(ν) dν = f(λ) dλ , dν = 1/λ , f(ν) 1/λ^2 dλ = f(λ) dλ, f(ν) = λ^2 f(λ)
    # S λ^2 f(λ) = σ
    # κ = σ * N
    # σ21 = S21 * λ210^2 # [m^3]
$γ = \left(\dfrac{T_{ref}}{T}\right)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})

Niso = N * NCO2 * iso_c

λ_ul = λ_ul0 / (1.0 + λ_ul0 * δ_a * p)

# γ = (Tref/T)^n_{air} (γ_a(p_{ref], T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref], T_{ref})
dT = (TREF/T)^n_a
γ = dT * (γ_a * p * (1.0 - NCO2) + γ_s * p * NCO2)

ΔλL = λ_ul^2 * γ
ΔλG = sqrt(2.0 * c_kB * T / iso_m) / c_c * λ_ul
ΔλL_mean[iλ] = ΔλL
ΔλD_mean[iλ] = ΔλG

$N_l  = \dfrac{g_l}{Q(T, iso)} \exp(- E_l  β)  N_{iso}$
$N_u  = \dfrac{g_u}{Q(T, iso)} \exp(- E_u  β)  N_{iso}$

$ϵ(λ) = \dfrac{h c}{λ_0} N_u A_{ul} * f(λ) \dfrac{dΩ}{4 π}$
$κ(λ) = \dfrac{h c}{λ_0} N_l B_{lu}  \left(1 - \dfrac{N_u}{N_l}  \dfrac{g_l}{g_u}\right)  \dfrac{λ_0^2}{c} f(λ)$
