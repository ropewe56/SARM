module Sarm

include("input/atmosphere.jl")
include("input/moleculardata.jl")
include("input/paths.jl")
include("input/runparameter.jl")

include("utils.jl")
include("planck.jl")
include("profiles.jl")
include("profiles2.jl")

include("resultdata.jl")
include("run_spectrum.jl")
include("spectrum_v1.jl")
include("start_julia.sh")

include("spectrum.jl")
include("main.jl")

end