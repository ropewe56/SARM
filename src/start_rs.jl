# starting via shared lib
#   CURL_SSL_VERSION mismatch 
#   mv  /home/wester/.julia/juliaup/julia-1.11.6+0.x64.linux.gnu/bin/../lib/julia/libcurl.so.4.8.0 \
#       /home/wester/.julia/juliaup/julia-1.11.6+0.x64.linux.gnu/bin/../lib/julia/libcurl.so.4.8.0_backup
#   ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4.8.0 \
#         /home/wester/.julia/juliaup/julia-1.11.5+0.x64.linux.gnu/lib/julia 

const SAMR_SO = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/target/release/libsarm.so"

function start_rust(subdir; so=false)
    jsonfile = joinpath(RESULT_RS, subdir, "parameter.json")
     
    if so
        @ccall SAMR_SO.init_jl(jsonfile::Cstring)::Cvoid
    else
        sarm_bin = "/home/wester/Projects/Julia/Climate-Energy/Sarm.rs/target/release/sarm"
        run(`$sarm_bin  $jsonfile`)
    end
end
