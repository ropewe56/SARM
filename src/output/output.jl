using Printf
using SpecialFileIO

"""
    Save convolved absorption and emission spectra to hdf5 file
"""
function save_convolved_to_hdf5(λb, κb, ϵb, hdf5_path)
    groups = Dict("wl_bin" => Dict("wl" => λb, "kappa" => κb, "epsilon" => ϵb))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function save_spectrum_as_hdf5(hdf5_path1, hdf5_path2, λb, Iλb, κb, ϵb, md)
    path1 = @sprintf("%s%s", plitext(hdf5_path)[1], ";spec.hdf5")
    path2 = @sprintf("%s%s", plitext(hdf5_path)[1], ";atm.hdf5")

    groups_spec = Dict("spec" => Dict("wl" => λb, "kappa" => κb, "eps" => ϵb))
    save_groups_as_hdf5(path1, groups_spec; permute_dims_p=false, extension=".hdf5", script_dir=false)    

end

function write_to_hdf5(paths, atm, outid, ih, iN, iθ, md, λb, Iλb, κb, ϵb)
    atmname = @sprintf("atm_%03d_%03d_%d_%d_%5.3e", outid, ih, iN, iθ, atm.h[ih])
    hdf5_path = joinpath(paths.outdir, "atm", atmname)

    groups_atm = Dict("groups_atm" => Dict("h" => atm.h, "p" => atm.p, "T" => atm.T, "N" => atm.N))
    save_groups_as_hdf5(hdf5_path, groups_atm; permute_dims_p=false, extension=".hdf5", script_dir=false)    

    sepctrumname = @sprintf("spectrum_%03d_%03d_%d_%d_%5.3e", outid, ih, iN, iθ, atm.h[ih])
    hdf5_path = joinpath(paths.outdir, "spectrum", sepctrumname)
    groups_spec = Dict( "Ib"         => Dict("wl"    => λb,     "Ib" => Iλb), 
                        "kappa_eps1" => Dict("kappa" => κb[1], "eps" => ϵb[1]),
                        "kappa_eps2" => Dict("kappa" => κb[2], "eps" => ϵb[2]))
        
    save_groups_as_hdf5(hdf5_path, groups_spec; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end