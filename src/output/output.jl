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
    hdf5_path = joinpath(paths[:atm], "atm.hdf5")
    groups = Dict("atm" => Dict("h" => atm.h, "p" => atm.p, "T" => atm.T, "N" => atm.N))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function write_results_to_hdf5(paths, atm, ic, iθ, ih, ML, λb, Iλb, κb, ϵb, κbs, ϵbs)
    spectrum_name = @sprintf("spectrum_%03d_%d_%d_%4.1f.hdf5", ic, iθ, ih, atm.h[ih]*1.0e-3)
    hdf5_path = joinpath(paths[:spectrum], spectrum_name)

    # ML :     iso, S21, λ21, γ, ΔλL, ΔλG, N1, N2, miso[iso], ϵ, κ1, κ2

    d = Dict("λ" => λb, "I" => Iλb, "κ" => κb, "ϵ" => ϵb)
    for (spec, val) in κbs
        ks = @sprintf("κ_%s", spec)
        es = @sprintf("ϵ_%s", spec)
        d[ks] = κbs[spec]
        d[es] = ϵbs[spec]

        sl = @sprintf("Sl_%s", spec)
        ll = @sprintf("λl_%s", spec)
        el = @sprintf("ϵl_%s", spec)
        kl1 = @sprintf("κl1_%s", spec)
        kl2 = @sprintf("κl2_%s", spec)
        d[sl] = ML[spec][2,:] 
        d[ll] = ML[spec][3,:] 
        d[el] = ML[spec][10,:] 
        d[kl1] = ML[spec][11,:] 
        d[kl2] = ML[spec][12,:] 
    end

    groups = Dict( "sarm" => d)
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)
    hdf5_path
end