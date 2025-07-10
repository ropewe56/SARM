using Parameters
using PhysConst

include("../utils.jl")

struct Atmosphere
    h  :: Vector{Float64}
    p  :: Vector{Float64}
    T  :: Vector{Float64}
    N  :: Vector{Float64}
    ip :: Vector{Interpolations.Extrapolation}
end

function get_densities(atm, z, c0)
    c1 = H2O_concentration(z; c0=-1.0)
    c2 = CO2_concentration(z; c0=c0)
    N*c1, N*c2
end

#function CO2_concentration(z; c0 = 425.0)
#    h = [0.0, 10000.0, 70000.0]

"""
    Create a Vector{Float64} with z-values for the integration

    not used anymore, instead see make_z_log10
    zmin
    zmax
    dzmin
    dzmax
    n
    e
"""
function make_z_e(par)
    nz = par.nz
    e  = par.e
    dzmin, dzmax, zmin, zmax = par.dzmin, par.dzmax, par.zmin, par.zmax

    dz = collect(range(dzmin, dzmax, nz))
    dz = dz.^e
    z  = cumsum(dz)
    z  = z * zmax / maximum(z)

    zout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for zo in zout
        zz = z .- zo
        i = argmin(zz.^2)
        z[i] = zo
    end
    z 
end

function make_z_exp(par)
    nz = par.nz
    e  = par.e
    dzmin, dzmax, zmin, zmax = par.dzmin, par.dzmax, par.zmin, par.zmax

    dz = collect(range(dzmin, dzmax, nz))
    dz = dz.^e
    z  = cumsum(dz)
    z  = z * zmax / maximum(z)

    zout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for zo in zout
        zz = z .- zo
        i = argmin(zz.^2)
        z[i] = zo
    end

    h  = collect(range(zmin, zmax, par.nz))
    z = @. exp(h/maximum(h)*e) - 1.0
    z = z/maximum(z) * par.zmax
    z
end

"""
    make_z_log10(zmin, zmax, n)

    make a Vector of z-values: radiation transfer from z[i] to z[i+1]
    z-values are created equadistantly on a log10 scale, than z = 10^log10_z

    create a Vector{Int64} z_iout of length n, plot inetensity and spectrum where z_iou == 1

    zmin : starting z value (e.g. 0.1 m)
    zmax : end z value (70 km TAO)
    n : number of z-values

"""
function make_z_log10(par)
    log10_z = collect(range(log10(max(1.0e-1, par.xmin)), log10(par.zmax), par.nz))
    z = @. 10^log10_z
    z
end

function H2O_concentration(z; c0=-1.0)
    h = reverse([84.977, 76.278, 67.577, 32.608, 41.176, 52.132, 13.792, 11.565, 8.095, 6.142, 3.77, 1.952, 0.137].*1.0e3)
    c = reverse([-5.8683, -5.5885, -5.4074, -5.3251, -5.3086, -5.2757, -5.1934, -4.6173, -3.465, -3.0864, -2.642, -2.3951, -2.0988])
    index = sortperm(h)
    h2  = h[index]
    c2 = c[index]
    linear_interpolation(h2, c2, extrapolation_bc = Line())
end

function CO2_concentration(z; c0 = 425.0)
    h = [0.0, 10000.0, 70000.0]
    c = [1.0, 2.0/3.0, 2.0/3.0]
    linear_interpolation(h, c, extrapolation_bc = Line())
end

function get_concentration(atm, z, i, c0)
    c = 10.0.^(atm.ip[i](z))
    c0 = if c0 < 0.0
        c[1]
    end
    c ./ maximum(c) .* c0
end

"""
    get_zTpN(nz)

height dependent digitized values of T and p http://climatemodels.uchicago.edu/modtran/
    
create linear interpolation objects and interpolate values onto an equidistamt grid of n points

compute density N and inverse interpolate N,h with N equidistant

the h values are z values
"""
function get_pT_interpolator()
    # p over h
    h1 = [   0.0,   1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0,  11.0,  12.0,
            13.0,  14.0,  15.0,  16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 70.0]
    p1 = [1013.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 324.0, 281.0, 243.0, 209.0,
           179.0, 153.0, 130.0, 111.0, 95.0, 81.2, 69.5, 59.5, 51.0, 43.4, 27.7, 13.2, 6.52, 3.33, 1.76, 0.951, 0.067]
    x0  = h1 .* 1.0e3
    y0  = p1 .* 1.0e2
    
    ip_hp = linear_interpolation(x0, y0, extrapolation_bc = Line())

    # T over h
    x0 = [0.0, 13.0, 17.0, 25.0, 30.0, 45.0, 50.0, 70.0] .* 1.0e3
    y0 = [288.0, 215.8, 215.7, 225.1, 233.7, 269.9, 275.7, 218.1]
    ip_hT = linear_interpolation(x0, y0, extrapolation_bc = Line())

    ip_hp, ip_hT
end

function Atmosphere(par)
    ip_hp, ip_hT = get_pT_interpolator()

    z = if par.e == :e
        make_z_e(par)
    elseif par.zmethod == :exp
        make_z_exp(par)
    elseif par.zmethod == :log10
        make_z_log10(par)
    elseif par.zmethod == :equalnumber
        h = collect(range(par.zmin, par.zmax, par.nz*10))
        p, T, = ip_hp.(h), ip_hT.(h)
        N = @. p / (c_kB * T)

        np = par.nz*10
        x0 = sqrt.(reverse(N))
        y0 = reverse(h)
        x, y, ip_Nh = lininterp(x0, y0, np)
        z = reverse(y)
        N = reverse(x).^2
        z
    end

    p, T, = ip_hp.(z), ip_hT.(z)
    N = @. p / (c_kB * T)

    ip1 = H2O_concentration(z)
    ip2 = CO2_concentration(z)

    Atmosphere(z, p, T, N, [ip1, ip2])
end

#par = parameter()
#make_zpTN(par)