# CURL_SSL_VERSION
# (base) wester@HPLinux:~/.julia/juliaup/julia-1.11.5+0.x64.linux.gnu/lib/julia$ cp /lib/x86_64-linux-gnu/libcurl.so.4.8.0 .

jsonfile = "/home/wester/Projects/Julia/Climate-Energy/Sarm.jl/results/rs_CO2/parameter.json"
rslib = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/target/release/libsarm.so"
@ccall rslib.init_jl(jsonfile::Cstring)::Cvoid

sarm = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/target/release/sarm"
run(`$sarm  $jsonfile`)