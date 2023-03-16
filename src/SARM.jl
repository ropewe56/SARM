module SimpleRadTrans

include("create_input/create_QTp_lines_files.jl")
include("create_input/create_Nthetaz_json_files.jl")

include("parameter.jl")
include("profiles.jl")
include("resultdata.jl")
include("spectralline.jl")
include("planck.jl")

include("interpolator.jl")
include("spectrum.jl")

end