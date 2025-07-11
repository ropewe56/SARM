using Parameters
using PhysConst

include("../utils.jl")

struct Atmosphere
    h     :: Vector{Float64}
    p     :: Vector{Float64}
    T     :: Vector{Float64}
    N     :: Vector{Float64}
    chitp :: Vector{Interpolations.Extrapolation}
end

function get_densities(atm, h, c0)
    c1 = H2O_concentration(h; c0=-1.0)
    c2 = CO2_concentration(h; c0=c0)
    N*c1, N*c2
end

"""
    Create a Vector{Float64} with h-values for the integration

    not used anymore, instead see make_h_log10
    hmin
    hmax
    dhmin
    dhmax
    n
    e
"""
function make_h_e(par)
    nh = par.nh
    e  = par.e
    dhmin, dhmax, hmin, hmax = par.dhmin, par.dhmax, par.hmin, par.hmax

    dh = collect(range(dhmin, dhmax, nh))
    dh = dh.^e
    h  = cumsum(dh)
    h  = h * hmax / maximum(h)

    hout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for ho in hout
        hh = h .- ho
        i = argmin(hh.^2)
        h[i] = ho
    end
    h 
end

function make_h_exp(par)
    nh = par.nh
    e  = par.e
    dhmin, dhmax, hmin, hmax = par.dhmin, par.dhmax, par.hmin, par.hmax

    dh = collect(range(dhmin, dhmax, nh))
    dh = dh.^e
    h  = cumsum(dh)
    h  = h * hmax / maximum(h)

    hout = [0.1, 0.5, 1.0, 10.0, 100.0, 200.0, 500.0, 1000.0, 5000.0, 10000.0, 70000.0]
    for ho in hout
        hh = h .- ho
        i = argmin(hh.^2)
        h[i] = ho
    end

    h  = collect(range(hmin, hmax, par.nh))
    h = @. exp(h/maximum(h)*e) - 1.0
    h = h/maximum(h) * par.hmax
    h
end

"""
    make_h_log10(hmin, hmax, n)

    make a Vector of h-values: radiation transfer from h[i] to h[i+1]
    h-values are created equadistantly on a log10 scale, than h = 10^log10_h

    create a Vector{Int64} h_iout of length n, plot inetensity and spectrum where h_iou == 1

    hmin : starting h value (e.g. 0.1 m)
    hmax : end h value (70 km TAO)
    n : number of h-values

"""
function make_h_log10(par)
    log10_h = collect(range(log10(max(1.0e-1, par.xmin)), log10(par.hmax), par.nh))
    h = @. 10^log10_h
    h
end

function H2O_concentration(h)
    h = reverse([84.977, 76.278, 67.577, 32.608, 41.176, 52.132, 13.792, 11.565, 8.095, 6.142, 3.77, 1.952, 0.137].*1.0e3)
    c_log10 = reverse([-5.8683, -5.5885, -5.4074, -5.3251, -5.3086, -5.2757, -5.1934, -4.6173, -3.465, -3.0864, -2.642, -2.3951, -2.0988])

    index = sortperm(h)
    h2  = h[index]
    c_log10 = c_log10[index]

    c   = 10.0.^c_log10
    c1  = c ./ maximum(c)
    c2  = log10.(c1)
    linear_interpolation(h2, c2, extrapolation_bc = Line())
end

function CO2_concentration(h)
    h = [0.0, 10000.0, 70000.0]
    c = [1.0, 2.0/3.0, 2.0/3.0]
    linear_interpolation(h, c, extrapolation_bc = Line())
end

function get_concentrations(atm, h, ch0)
    ch1 = 10.0.^(atm.chitp[1](h))
    ch2 = 10.0.^(atm.chitp[2](h))

    ch1 = ch1 * ch0[1]
    ch2 = ch2 * ch0[2]

    [ch1, ch2]
end

"""
    get_hTpN(nh)

height dependent digitihed values of T and p http://climatemodels.uchicago.edu/modtran/
    
create linear interpolation objects and interpolate values onto an equidistamt grid of n points

compute density N and inverse interpolate N,h with N equidistant

the h values are h values
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

    h = if par.e == :e
        make_h_e(par)
    elseif par.hmethod == :exp
        make_h_exp(par)
    elseif par.hmethod == :log10
        make_h_log10(par)
    elseif par.hmethod == :equalnumber
        h = collect(range(par.hmin, par.hmax, par.nh*10))
        p, T, = ip_hp.(h), ip_hT.(h)
        N = @. p / (c_kB * T)

        np = par.nh*10
        x0 = sqrt.(reverse(N))
        y0 = reverse(h)
        x, y, ip_Nh = lininterp(x0, y0, np)
        h = reverse(y)
        N = reverse(x).^2
        h
    end

    p, T, = ip_hp.(h), ip_hT.(h)
    N = @. p / (c_kB * T)

    chitp1 = H2O_concentration(h)
    chitp2 = CO2_concentration(h)

    Atmosphere(h, p, T, N, [chitp1, chitp2])
end

#par = parameter()
#make_hpTN(par)