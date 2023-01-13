## $CO_2$ partial density
$N_{iso} = N(h) \;C_{CO2} \; C_{iso}$

## Pressure line shift
$λ_{ul} = \dfrac{λ_{ul0}}{1 + λ_{ul0}  \delta_{air}  p}$

## Line shape
$\int f(λ) d\lambda = 1$
### Lorentz line shape
$f_L(\lambda) = \dfrac{1}{\pi} \dfrac{\Delta \lambda_L^2}{(\lambda - \lambda_{ul})^2 + \Delta \lambda_L^2}$
$Δλ_L = γ λ_{ul}^2$
$γ = \left(\dfrac{T_{ref}}{T}\right)^{n_{air}} \left[γ_a(p_{ref}, T_{ref}) (p - p_{CO2}) + γ_s p_{CO2}(p_{ref}, T_{ref})\right]$
### Gauss line shape
$f_G(\lambda) = \sqrt{\dfrac{\ln 2}{\pi \Delta\lambda_G^2}} \exp\left(-\dfrac{\ln 2}{\Delta\lambda_G^2}(\lambda -\lambda_{ul})^2\right)$
$Δλ_G = \dfrac{λ_{ul}}{c} \sqrt{\dfrac{2  k_B  T}{M(iso)}}$

## Occupation densities
$N_l  = \dfrac{g_l}{Q(T, iso)} \exp\left(- \dfrac{E_l}{K_B T}\right)  N_{iso}$
$N_u  = \dfrac{g_u}{Q(T, iso)} \exp\left(- \dfrac{E_u}{K_B T}\right)  N_{iso}$

## Radiation transfer
$\dfrac{d^2 I(\lambda)}{ds \; d\lambda}  = -\dfrac{d I(\lambda)}{d\lambda}κ(λ) + \dfrac{d ϵ(λ)}{d\lambda}$

$\int \dfrac{d^2 I(\lambda)}{ds \; d\lambda} d\lambda = - \int \dfrac{d I(\lambda)}{d\lambda}κ(λ) d\lambda  + \int \dfrac{d ϵ(λ)}{d\lambda} d\lambda$

$\dfrac{d }{ds} \int \dfrac{d I(\lambda)}{d\lambda} d\lambda = - \int \dfrac{d I(\lambda)}{d\lambda}κ(λ) d\lambda  + \int \dfrac{d ϵ(λ)}{d\lambda} d\lambda$

$κ(λ) = \sum_i κ_i(λ)$

$\dfrac{dϵ(λ)}{d \lambda} = \sum_i \dfrac{dϵ_i(λ)}{d \lambda}$

$\int \dfrac{dϵ(λ)}{d \lambda} d \lambda = \sum_i \int \dfrac{dϵ_i(λ)}{d \lambda} d \lambda$

### Emission
$ϵ_i(λ) = \dfrac{h c}{λ_{ul}} N_u A_{ul} f_i(λ) \dfrac{dΩ}{4 π}$

### Absorption
$κ_i(λ) = \dfrac{h λ_{ul}}{c} \left(N_l B_{lu} - N_u B_{ul}\right) f_i(λ)$

$B_{ul} = \dfrac{1}{8 \pi} \dfrac{\lambda_{ul}^3}{h} SA_{ul}$

$B_{lu} = \dfrac{g_u}{g_l} B_{ul}$


## Intensity

$\epsilon(\lambda) = \sum_j \epsilon_j(\lambda)$
$\kappa(\lambda) = \sum_j \kappa_j(\lambda)$
$I(\lambda,s) = I(\lambda,s_0) \exp \left(- \kappa(\lambda) (s-s_0)\right) + \dfrac{\epsilon(\lambda)}{\kappa(\lambda)} \left(1 - \exp\left(-\kappa(\lambda) (s-s_0)\right)\right)$


## Line strength
* https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy

$A_{ki} = \dfrac{2 \pi e^2}{m_e c \epsilon_0 \lambda^2} \dfrac{g_i}{g_k} f_{ik}
= \dfrac{16  \pi^3}{3 h \epsilon_0 \lambda^3 g_k} S$
$\epsilon_{line} =  \dfrac{1}{4 \pi} h ν A_{ki} N_k  \left[\dfrac{J}{s m^3}\right]$

## Modran
   * http://climatemodels.uchicago.edu/modtran/
   * Spectrum Upward IR Flux
   * density
   * pressure
   * CO2 concentration
