using PyPlot
pygui(false)
pygui(:qt5)

using Common
using PhysConst

function plot_lines()
    input_file = joinpath(@__DIR__, "CO2.hdf5")
    group_name, dataset_name = "CO2", "data"
    data = read_npy(input_file)

    # 0 wlnm
    # 1 dEnm
    # 2 Em
    # 3 Anm
    # 4 gn
    # 5 gm

    n1, n2 = size(data)

    wl = data[1,:]
    Em = data[2,:]
    En = Em + data[3,:]
    A  = data[4,:]
    gn = data[5,:]
    L  = data[12,:]

    T  = 288.0
    c1 = 1.0/(4.0*π) * c_h * c_c

    plot(wl*1.0e6, np.exp(-En / (c_kB * T)) * A * gn * c1 / 2.900e2 /wl * 2.5e23 * 400.0e-6 * 1.0e3)
    figure()

    input_file = joinpath(@__DIR__, "ea.hdf5")
    group = "wea"
    dataset = "wea"
    wea = load_array_as_hdf5(input_file; script_dir=false)

    plot(wea[:,0]*1.0e6,wea[:,1])
    figure()
    plot(wea[:,0]*1.0e6,wea[:,2])
end

function plot_intensity_vs_wavelength(ifig, path, i)
    lI  = load_array_as_hdf5(joinpath(path[1],  string(path[2], path[3])), group="intensity", dataset="l_I")

    @infoe size(lI)
    figure()
    plot(lI[1,:]*1.0e6, lI[2,:])
    xlabel("λ [μm]")
    ylabel("I_λ [W/m^2/sr/m]")
    savefig(joinpath(path[1], "png", "I", @sprintf("%d_%s_I_λ.png", i, path[2])))
    close(ifig)
end

function compose_path(path)
    joinpath(path[1], string(path[2], path[3]))
end

function plot_intensity_over_wavelength(ifig, planck_path, planck_intensity_part, path, int_nb, istart)
    group_name, dataset_name = "A", "A"

    lIPlanck = load_array_as_hdf5(compose_path(planck_path), group="A", dataset="A")
    lIPlanckPart = load_array_as_hdf5(compose_path(planck_intensity_part), group="A", dataset="A")

    lI = load_array_as_hdf5(compose_path(path), group="A", dataset="A")

    figure(ifig)
    plot(lI[1,:]*1.0e6, lI[2,:], "b")
    xlabel("λ [μm]")
    ylabel("I_λ [W/m^2/sr/m]")
    pngdir = joinpath(path[1], "png", "I", "noplanck")
    mkpath(pngdir)
    png_path = joinpath(pngdir, @sprintf("%03d_%s_I_wl.png", int_nb, path[2]))
    savefig(png_path)
    close(ifig)

    figure(ifig+1)
    plot(lI[1,:]*1.0e6, lI[2,:], "b")
    plot(lIPlanck[1,:]*1.0e6, lIPlanck[2,:], "r")
    xlabel("λ [μm]")
    ylabel("I_λ [W/m^2/sr/m]")
    pngdir = joinpath(path[1], "png", "I", "planck_1")
    mkpath(pngdir)
    png_path = joinpath(pngdir, @sprintf("%03d_%s_I_wl.png", int_nb, path[2]))
    savefig(png_path)
    close(ifig+1)

    figure(ifig+2)
    plot(lI[1,:]*1.0e6, lI[2,:], "b")
    n1, n2 = size(lIPlanckPart)
    for i in 2:n1
      #@infoe i, n1, n2
      plot(lIPlanckPart[1,:]*1.0e6, lIPlanckPart[i,:])
    end
    xlabel("λ [μm]")
    ylabel("I_λ [W/m^2/sr/m]")
    pngdir = joinpath(path[1], "png", "I", "planck_2")
    mkpath(pngdir)
    png_path = joinpath(pngdir, @sprintf("%03d_%s_I_wl.png", int_nb, path[2]))
    savefig(png_path)

    close(ifig+2)
    istart
end

function plot_spectrum_over_wavelength(ifig1, ifig2, path, i)
    group_name, dataset_name = "conv", "conv"
    hdf5path = joinpath(path[1], string(path[2], path[3]))
    @infoe hdf5path
    wl = load_array_as_hdf5(hdf5path, group="A", dataset="A")

    λ = wl[1,:]*1.0e6  # µm
    κ = wl[2,:]        # 1/m
    ϵ = wl[3,:]        # W / m^3 / sr / m

    pngdir = joinpath(path[1], "png", "S")
    mkpath(pngdir)
    pngpath = joinpath(pngdir, @sprintf("%03d_%s_conv_kap.png", i, path[2]))
    #@infoe @sprintf("%s  %s", size(wl), pngpath)

    figure(ifig1)
    plot(λ, κ, "b")
    xlabel("λ [μm]")
    ylabel("κ [1/m]")
    savefig(pngpath)
    close(ifig1)

    figure(ifig2)
    plot(λ, ϵ, "b")
    xlabel("λ [μm]")
    ylabel("ϵ [W/m^3/sr/m]")
    savefig(joinpath(pngdir, @sprintf("%03d_%s_conv_eps.png", i, path[2])))
    close(ifig2)
end

function renames(root, ext)
    spectrum, intensity, moving_average_k, planck_intensity = [], [], nothing, nothing
    dirs = Set()
    for (rr, dd, ff) in walkdir(root)
        for f in ff
            fs = splitext(f)
            if fs[2] == ext
                if "ma_k" in f
                    a = "ma_k"
                    b = "moving_average_k"
                    shutil.move(joinpath(rr, fs[1]+fs[2]), joinpath(rr, fs[1].replace(a, b)+fs[2]))
                elseif "lI_planck" in f
                    a = "lI_planck"
                    b = "planck_intensity"
                    shutil.move(joinpath(rr, fs[1]+fs[2]), joinpath(rr, fs[1].replace(a, b)+fs[2]))
                elseif "conv" in f
                    a = "conv"
                    b = "spectrum"
                    shutil.move(joinpath(rr, fs[1]+fs[2]), joinpath(rr, fs[1].replace(a, b)+fs[2]))
                elseif "lI" in f
                    a = "lI"
                    b = "intensity"
                    shutil.move(joinpath(rr, fs[1]+fs[2]), joinpath(rr, fs[1].replace(a, b)+fs[2]))
                end
            end
        end
    end
end

function get_all_files_by_extension(root, ext)
    ispectrum, iintensity = [], []
    spectrum, intensity, moving_average_k, planck_intensity, planck_intensity_part = [], [], nothing, nothing, nothing
    dirs = Set()
    for (rr, dd, ff) in walkdir(root)
        for f in ff
            try
                fs = splitext(f)
                #@infoe @sprintf("%s | %s | %s", f, fs, ext)
                if fs[2] == ext
                    if occursin("moving_average_k", f)
                        moving_average_k = (rr, fs[1], fs[2])
                    elseif occursin("planck_intensity_part", f)
                        planck_intensity_part = (rr, fs[1], fs[2])
                    elseif occursin("planck_intensity", f)
                        planck_intensity = (rr, fs[1], fs[2])
                    elseif occursin("result_", f)
                        zs = split(fs[1], "_")
                        #@infoe zs
                        iN = floor(Int64, parse(Float64, zs[2]))
                        iθ = floor(Int64, parse(Float64, zs[3]))
                    else
                        zs = split(fs[1], "_")
                        iN = floor(Int64, parse(Float64, zs[2]))
                        iθ = floor(Int64, parse(Float64, zs[3]))
                        iz = floor(Int64, parse(Float64, zs[4]))
                        #@infoe @sprintf("%8d  %s %s %s", iz, rr, fs[1], fs[2])
                        if occursin("spectrum", f)
                            push!(ispectrum, (iN, iθ, iz))
                            push!(spectrum, (rr, fs[1], fs[2]))
                            push!(dirs, rr)
                        elseif occursin("intensity", f)
                            push!(iintensity, (iN, iθ, iz))
                            push!(intensity, (rr, fs[1], fs[2]))
                        end
                    end
                end
            catch e
                @warne e
                @warne f
            end
        end
    end

    index = sortperm(ispectrum)
    spectrum = spectrum[index]
    index = sortperm(iintensity)
    intensity = intensity[index]

    spectrum, intensity, moving_average_k, planck_intensity, planck_intensity_part
end

function plot_moving_average(path)
    ma = read_npy(joinpath(path[1], path[2]+path[3]))
    plot(ma[:,1]*1.0e6, ma[:,2])
    figure()
end


function go(subdir)
    out_root = joinpath(dirname(@__DIR__), "radoutput")

    subsubdirs = readdir(joinpath(out_root, subdir))
    sort!(subsubdirs)

    #for d in subsubdirs:
    #    renames(joinpath(out_root, subdir, d), ".npy")

    istart = 0
    for d in subsubdirs
        spectrum, intensity, moving_average_k, planck_intensity, planck_intensity_part = get_all_files_by_extension(joinpath(out_root, subdir, d), ".hdf5")
        @infoe planck_intensity_part
        if !(intensity === nothing) && !(planck_intensity === nothing)
            n = size(intensity,1)
            for int_nb in 1:n
                istart = plot_intensity_over_wavelength(3, planck_intensity, planck_intensity_part, intensity[int_nb], int_nb, istart)
                plot_spectrum_over_wavelength(1, 2, spectrum[int_nb], int_nb)
            end
        end
    end
end
subdir = "A"
go(subdir)

