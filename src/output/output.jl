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

function write_atmosphere_to_hdf5(paths, atm)
    hdf5_path = joinpath(paths[:atm], "atm.hdf5")
    groups = Dict("atm" => Dict("h" => atm.h, "p" => atm.p, "T" => atm.T, "N" => atm.N))
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)    
end

function write_results_to_hdf5(paths, atm, ic, iθ, ih, linedata_dict, λb, Iλb, κb, ϵb, κbs, ϵbs, intfs)
    spectrum_name = @sprintf("spectrum_%d_%d_%03d_%07.1f.hdf5", ic, iθ, ih, atm.h[ih])
    hdf5_path = joinpath(paths.spectrum, spectrum_name)

    # linedata
    # 1             2          3          4          5    6    7   8    9    10  11  12 13
    # Float64(iso), miso[iso], aiso[iso], Qiso[iso], S21, λ21, γp, ΔλL, ΔλG, N1, N2, ϵ, κ

    d = Dict("λb" => λb, "Iλb" => Iλb, "κb" => κb, "ϵb" => ϵb)
    for (spec, val) in κbs
        linedata = linedata_dict[spec]

        ks = @sprintf("κb_%s", spec)
        es = @sprintf("ϵb_%s", spec)
        fs = @sprintf("intf_%s", spec)
        d[ks] = κbs[spec]
        d[es] = ϵbs[spec]
        d[fs] = intfs[spec]

        sl     = @sprintf("Sl_%s", spec)
        ll     = @sprintf("λl_%s", spec)
        el     = @sprintf("ϵl_%s", spec)
        kl     = @sprintf("κl_%s", spec)
        d[sl]  = [linedata[i][ 5] for i in eachindex(linedata)]
        d[ll]  = [linedata[i][ 6] for i in eachindex(linedata)]
        d[el]  = [linedata[i][12] for i in eachindex(linedata)]
        d[kl]  = [linedata[i][13] for i in eachindex(linedata)]
    end

    groups = Dict( "sarm" => d)
    save_groups_as_hdf5(hdf5_path, groups; permute_dims_p=false, extension=".hdf5", script_dir=false)
    hdf5_path
end