using Interpolations

const PPM = 1.0e-6
const C0H2O_PPM = 7966.0
const MOLECULE_SYMBOLS = [:H2O, :CO2]

"""
    CH2O_concentration(hi)
    normalized concentraion at height points hi
"""
function H2O_normalized_concentration_over_h(hi)
    # digitized values
    h = reverse([84.977, 76.278, 67.577, 32.608, 41.176, 52.132, 13.792, 11.565, 8.095, 6.142, 3.77, 1.952, 0.137].*1.0e3)
    c_log10 = reverse([-5.8683, -5.5885, -5.4074, -5.3251, -5.3086, -5.2757, -5.1934, -4.6173, -3.465, -3.0864, -2.642, -2.3951, -2.0988])

    index = sortperm(h)
    h2    = h[index]
    cl = c_log10[index]
    itp = linear_interpolation(h2, cl, extrapolation_bc = Line())

    cli = itp(atm.h)
    ci = 10.0.^cli
    ci = ci ./ maximum(ci)
    ci
end

"""
    CO2_concentration(hi)
    normalized concentraion at height points hi
"""
function CO2_normalized_concentration_over_h(hi)
    h = [0.0, 10000.0, 70000.0]
    c = [1.0, 2.0/3.0, 2.0/3.0]
    itp = linear_interpolation(h, c, extrapolation_bc = Line())
    itp(hi)
end

"""
    get_normalized_molcule_concentration_over_h(molecule, hi)
"""
function get_normalized_molecule_concentration_over_h(species, hi)
    ci = if species == :H2O
        H2O_normalized_concentration_over_h(hi)
    elseif species == :CO2
        CO2_normalized_concentration_over_h(hi)
    else
        @error species, "not implemented"
    end
    ci
end

"""
    get_concentrations(moleculardata, ih)
"""
function get_concentrations_over_h(moleculardata::Dict{Symbol,MolecularData}, ch0::Dict{Symbol, Vector{Float64}})
    ch = Dict{Symbol, Matrix{Float64}}()
    nbc = 0
    for (species, mddata) in moleculardata
        cnh  = mddata.cnh # Vector{Float}
        nbc = max(nbc, length(ch0[species]))
        M = cnh * (ch0[species] .* PPM)'
        ch[species] = M
    end
    ch, nbc
end

