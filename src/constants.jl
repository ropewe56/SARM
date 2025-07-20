using PhysConst

const hc = c_h * c_c

const LOG2 = log(2.0)

const TREF = 296.0 # https://hitran.org/docs/definitions-and-units/

const PPM = 1.0e-6
const C0H2O_PPM = 7966.0
const MOLECULE_SYMBOLS = [:H2O, :CO2]

const OUTROOT = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results"

const _cm = 1.0e-2
const _atm = 1.01325e5