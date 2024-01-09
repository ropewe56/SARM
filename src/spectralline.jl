using StaticArrays
using StatsBase
using HDF5

using RWPhysConst

const i_λ_ul0 = 1
const i_E_l   = 2
const i_E_u   = 3
const i_S     = 4
const i_A_ul  = 5
const i_γ_a   = 6
const i_γ_s   = 7
const i_n_a   = 8
const i_δ_a   = 9
const i_g_u   = 10
const i_g_l   = 11
const i_B_ul  = 12
const i_B_lu  = 13
const i_ΔλL0  = 14
const i_iso_m = 15
const i_iso_c = 16

const i_γ   = 1
const i_ΔλL = 2
const i_ΔλG = 3
const i_N_l = 4
const i_N_u = 5
const i_ϵ   = 6
const i_κ   = 7
const i_λ_ul= 8

struct LineData
    iso_ids :: Vector{Int64}
    linesc  :: Matrix{Float64} ##Vector{SVector{16, Float64}}
    linesv  :: Matrix{Float64} ##Vector{SVector{8, Float64}}
    nb_lines:: Int64
end

function LineData(specdat_path, max_isotope_id)
    # spectral data (HITRAN)
    spectral_data = read_npy(specdat_path)
    @infoe typeof(spectral_data), size(spectral_data)

    n1, n2 = size(spectral_data)
    if n1 > n2
        spectral_data = colm_to_rowm(spectral_data)
    end

    nb_lines = size(spectral_data,2)

    @infoe @sprintf("max_isotope_id = %d", max_isotope_id)
    iso_ids = Vector{Int64}(undef,0)
    for iline in 1:nb_lines
        iso_id = spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id
            push!(iso_ids, isotope_id)
        end
    end
    @infoe @sprintf("min_iso_id = %d, max_iso_id = %d, size(iso_ids,1) = %d", minimum(iso_ids), maximum(iso_ids), size(iso_ids,1))

    linesc = Matrix{Float64}(undef, 16, size(iso_ids,1))
    i = 1
    for iline in 1:nb_lines
        iso_id = spectral_data[12, iline]
        isotope_id = floor(Int64, iso_id)
        if isotope_id <= max_isotope_id

            λ_ul0  = spectral_data[ 1, iline]
            E_l    = spectral_data[ 2, iline]
            E_u    = spectral_data[ 3, iline]
            S      = spectral_data[ 4, iline]
            A_ul   = spectral_data[ 5, iline]
            γ_a    = spectral_data[ 6, iline]
            γ_s    = spectral_data[ 7, iline]
            n_a    = spectral_data[ 8, iline]
            δ_a    = spectral_data[ 9, iline]
            g_u    = spectral_data[10, iline]
            g_l    = spectral_data[11, iline]
            iso_m  = spectral_data[13, iline]
            iso_c  = spectral_data[14, iline]

            # Einstein coefficient of induced emission
            B_ul = A_ul * λ_ul0^3 / (8.0*π * c_h)
            # Einstein coefficient of absorption
            B_lu = g_u / g_l * B_ul;

            # Line pressure broadening coeffcient
            # Lorentz line width
            ΔλL0 = λ_ul0^2 * (γ_a * 1.0e5)

            linesc[i_λ_ul0, i] = λ_ul0
            linesc[i_E_l  , i] = E_l
            linesc[i_E_u  , i] = E_u
            linesc[i_S    , i] = S
            linesc[i_A_ul , i] = A_ul
            linesc[i_γ_a  , i] = γ_a
            linesc[i_γ_s  , i] = γ_s
            linesc[i_n_a  , i] = n_a
            linesc[i_δ_a  , i] = δ_a
            linesc[i_g_u  , i] = g_u
            linesc[i_g_l  , i] = g_l
            linesc[i_B_ul , i] = B_ul
            linesc[i_B_lu , i] = B_lu
            linesc[i_ΔλL0 , i] = ΔλL0
            linesc[i_iso_m, i] = iso_m
            linesc[i_iso_c, i] = iso_c

            i += 1
        end
    end

    nb_lines = size(iso_ids,1)
    linesv = Matrix{Float64}(undef, 8, nb_lines)
    LineData(iso_ids, linesc, linesv, nb_lines)
end

function save_lines_data_as_hdf5(ld::LineData, fname)
    @infoe fname
    ldc = ld.linesc
    ldv = ld.linesv
    a = Matrix{Float64}(undef, 24, nb_lines)
    for iline in 1:nd.nb_lines
        a[ 1, iline] = ldc[i_λ_ul0]
        a[ 2, iline] = ldc[i_E_l  ]
        a[ 3, iline] = ldc[i_E_u  ]
        a[ 4, iline] = ldc[i_S    ]
        a[ 5, iline] = ldc[i_A_ul ]
        a[ 6, iline] = ldc[i_γ_a  ]
        a[ 7, iline] = ldc[i_γ_s  ]
        a[ 8, iline] = ldc[i_n_a  ]
        a[ 9, iline] = ldc[i_δ_a  ]
        a[10, iline] = ldc[i_g_u  ]
        a[11, iline] = ldc[i_g_l  ]
        a[12, iline] = ldc[i_B_ul ]
        a[13, iline] = ldc[i_B_lu ]
        a[14, iline] = ldc[i_ΔλL0 ]
        a[15, iline] = ldc[i_iso_m]
        a[16, iline] = ldc[i_iso_c]
        a[17, iline] = ldv[i_γ   ]
        a[19, iline] = ldv[i_ΔλL ]
        a[20, iline] = ldv[i_ΔλG ]
        a[21, iline] = ldv[i_N_l ]
        a[22, iline] = ldv[i_N_u ]
        a[23, iline] = ldv[i_ϵ   ]
        a[24, iline] = ldv[i_κ   ]
    end
    save_array_as_hdf5(a, fname, script_dir=false)
end

function wavelengths_histogram(ld::LineData, nbins, hdf5_path)
    @infoe hdf5_path
    ldc = ld.linesc

    λmin = minimum(ldc[i_λ_ul0,:])
    λmax = maximum(ldc[i_λ_ul0,:])
    Δλ = (λmax - λmin) / Float64(nbins)
    edges = λmin:Δλ:λmax

    iso_min = minimum(ld.iso_ids)
    iso_max = maximum(ld.iso_ids)
    @infoe @sprintf("iso_min = %d, iso_max = %d", iso_min, iso_max)

    local gid
    group = "wavelengths_histogram"
    h5open(hdf5_path, "w") do fid
        try
            gid = g_create(fid, group)
        catch
            gid = create_group(fid, group)
        end

        h = fit(Histogram, ldc[i_λ_ul0,:], edges)
        gid["edges"] = collect(h.edges[1])
        gid["weights"] = h.weights
        @infoe @sprintf("nblines = %d", sum(h.weights))

        for iso in iso_min:iso_max
            index = @. ifelse(ld.iso_ids == iso, true, false)
            λ_ul0 = ldc[i_λ_ul0,index]
            h = fit(Histogram, λ_ul0, edges)
            dataset = @sprintf("weights_%d", iso)
            @infoe @sprintf("iso = %d, %s : size = %s, nblines = %d", iso, dataset, size(h.weights), sum(h.weights))
            gid[dataset] = h.weights
        end
    end
end

function energy_energy_histogram(ld::LineData, lnbins, Enbins, hdf5_path)
    @infoe hdf5_path
    ldc = ld.linesc

    λmin = minimum(ldc[i_λ_ul0,:])
    λmax = maximum(ldc[i_λ_ul0,:])

    Δλ = (λmax - λmin) / Float64(lnbins)
    λedges = λmin:Δλ:λmax

    Emin = minimum(ldc[i_E_l,:])
    Emax = maximum(ldc[i_E_l,:])
    ΔE = (Emax -Emin) / Float64(Enbins)
    Eedges = Emin:ΔE:Emax

    iso_min = minimum(ld.iso_ids)
    iso_max = maximum(ld.iso_ids)

    @infoe @sprintf("ee-hist: λmin    = %8.2e, λmax    = %8.2e", λmin, λmax)
    @infoe @sprintf("ee-hist: Emin    = %8.2e, Emax    = %8.2e", Emin, Emax)
    @infoe @sprintf("ee-hist: iso_min = %d, iso_max = %d", iso_min, iso_max)

    h = fit(Histogram, (ldc[i_λ_ul0,:], ldc[i_E_l,:]), (λedges, Eedges))

    @infoe @sprintf("ee-hist: size(h.weights) = %s", size(h.weights))

    local gid
    group = "wavelengths_histogram"
    h5open(hdf5_path, "w") do fid
        try
            gid = g_create(fid, group)
        catch
            gid = create_group(fid, group)
        end

        gid["ledges"] = collect(h.edges[1])
        gid["Eedges"] = collect(h.edges[2])
        gid["weights"] = h.weights
    end
    h.edges, h.weights
end
