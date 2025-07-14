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

function write_atm_to_hdf5(paths, atm)
    hdf5_path = joinpath(paths.results, "atm.hdf5")
    groups = Dict("atm" => Dict("h" => atm.h, "p" => atm.p, "T" => atm.T, "N" => atm.N))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function write_results_to_hdf5(paths, atm, outid, ih, iN, iθ, λb, Iλb, κbs, ϵbs)
    spectrum_name = @sprintf("spectrum_%03d_%03d_%d_%d_%4.1f.hdf5", outid, ih, iN, iθ, atm.h[ih]*1.0e-3)
    hdf5_path = joinpath(paths.results, "spectrum", spectrum_name)

    d = Dict("wl" => λb, "I" => Iλb)
    for (spec, val) in κbs
        k = @sprintf("k_%s", spec)
        e = @sprintf("e_%s", spec)
        d[k] = κbs[spec]
        d[e] = ϵbs[spec]
    end
    groups = Dict( "I" => d)
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end