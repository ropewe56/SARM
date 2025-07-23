using PhysConst.UnitConst

w21   = uconvert(u"m^(-1)", 2285.117634/u_cm)
S21r  = uconvert(u"m", 4.63E-30 * u_cm)      
A21   = 393.2 / u_s                          
γair  = uconvert(u"m^(-1)*Pa^(-1)", 0.0702 / u_cm / (1.0e5*u_Pa))
γself = uconvert(u"m^(-1)*Pa^(-1)", 0.092 / u_cm / (1.0e5*u_Pa))
w1    = uconvert(u"m^(-1)", 2586.2758 / u_cm)
nair  = 0.7
δair  = uconvert(u"m^(-1)*Pa^(-1)", -0.003055/ u_cm / u_Pa)
g2    = 1071
g1    = 1113

hc    = u_h*u_c
λ1    = 1.0/w1
λ21   = 1.0/w21

T = 300.0 * u_K
p = 1.00e5 * u_Pa
TREF = 288.0 * u_K
β  = 1.0/(u_kB * T)
βr = 1.0/(u_kB * TREF)

ciso = 0.004
N = 1.0e22 / u_m^3
Niso = ciso * N                                                         #  [1/m^3]

# frequency
ν21   = u_c/λ21
ν1    = u_c/λ1
B21   = A21 * u_c^3 / (8π * u_h * ν21^3)
B12   = g2 / g1 * B21

λ21 = λ / (1.0 + δair * λ * p)


γp = (TREF/T)^nair * (γair * p * (1.0 - ciso) + γself * p * ciso)      # [1/m]
ΔλL = λ21^2 * γp                                                 # [m]


κ21   = u_h/u_c * B21 * (g2/g1*N1 - N2) * ν21 * 1.0e-4
ϵ21   = u_h * ν21 * A21 * N2 * fν 

fλ = 1.0 / (π * ΔλL * (1.0 + ((λ21 - λ0)/ΔλLh)^2))



# Doppler broadening
m = 1.0e-27*u_kg
a = uconvert(NoUnits, sqrt(2.0 * u_kB * T / (m * u_c^2)))
ΔλG = uconvert(u"m", a *  λ21)

# occupation numbers
E1 = u_h * u_c/λ1
E2 = u_h * u_c/λ1 + u_h * u_c/λ21
N1 = g1 * exp(- E1 * β) * Niso
N2 = g2 * exp(- E2 * β) * Niso

ϵ = uconvert(u"W * m^-(3)", hc/λ21 * N2 * A21 * 1.0 / (4.0 * π))

κ = u_h * λ21 / u_c * (N1 * B12 - N2 * B21)                      # Js * m * s/m / m^3 * m^3/(J*s^2) = 1

ν = u_c/λ21

F = 2*u_h * ν^3 / u_c^3

B21 = A21/F
B12 = B21

κ_ = u_h * u_c/λ21 * (N1 * B12 - N2 * B21) / ΔλL
uconvert(u"m^3/J/s^2", B12)

B12 = uconvert(u"m^3/J/s^2", u_e^2/(4u_ϵ0 * u_me * u_h * ν))
B21 = B12
uconvert(u"m^3/J/s^2", B12)

uconvert(NoUnits, 1.0u"μm/m")

ν = c/λ
dν = u_c/λ21^2 * dλ