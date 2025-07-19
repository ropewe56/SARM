using Parameters
using PhysConst

include("../utils.jl")

struct Atmosphere
    h_iout:: Vector{Int64}
    h     :: Vector{Float64}
    p     :: Vector{Float64}
    T     :: Vector{Float64}
    N     :: Vector{Float64}
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
    nh = par[:nh]
    e  = par[:e]
    dhmin, dhmax, hmin, hmax = par[:dhmin], par[:dhmax], par[:hmin], par[:hmax]

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
    nh = par[:nh]
    e  = par[:e]
    dhmin, dhmax, hmin, hmax = par[:dhmin], par[:dhmax], par[:hmin], par[:hmax]

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

    h  = collect(range(hmin, hmax, par[:nh]))
    h = @. exp(h/maximum(h)*e) - 1.0
    h = h/maximum(h) * par[:hmax]
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
    log10_h = collect(range(log10(max(1.0e-1, par[:xmin])), log10(par[:hmax]), par[:nh]))
    h = @. 10^log10_h
    h
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
    h1 = h1 .* 1.0e3
    p1 = p1 .* 1.0e2    
    ip_hp = linear_interpolation(h1, p1, extrapolation_bc = Line())

    # T over h
    h2 = [0.0, 13.0, 17.0, 25.0, 30.0, 45.0, 50.0, 70.0] .* 1.0e3
    T2 = [288.0, 215.8, 215.7, 225.1, 233.7, 269.9, 275.7, 218.1]
    ip_hT = linear_interpolation(h2, T2, extrapolation_bc = Line())

    ip_hp, ip_hT, h1, p1, h2, T2
end

function Atmosphere(par)
    ip_hp, ip_hT, h1, p1, h2, T2 = get_pT_interpolator()
    T12 = ip_hT.(h1)
    N12 = @. p1 / (c_kB * T12)

    h = if par[:e] == :e
        make_h_e(par)
    elseif par[:hmethod] == :exp
        make_h_exp(par)
    elseif par[:hmethod] == :log10
        make_h_log10(par)
    elseif par[:hmethod] == :equalnumber
        np = par[:nh]*100

        h = collect(range(par[:hmin], par[:hmax], np))
        p = ip_hp.(h)
        T = ip_hT.(h)
        N = @. p / (c_kB * T)

        # >> number of particles within a bin
        hh = @. 0.5 * (h[2:end] + h[1:end-1])
        dh = (h[2:end] - h[1:end-1])
        Nh = @. 0.5*(N[2:end] + N[1:end-1]) * dh
        h2, N2, ip = lininterp(hh, Nh, np)
        h3 = collect(range(h[1], h[end], np)) 
        N3 = ip.(h3)
        # << number of particles within a bin

        # >> interpolate N, h
        Ni, hi, ip = lininterp(reverse(N3), reverse(h3), np) # knot vectors must be uinique and increasing
        # equidistant number of particles within a bin
        x1, x2 = 0.05, 3.0
        ff = collect(range(x1, 1.0, par[:nh])).^x2
        fff = @. (ff - ff[1]) / (ff[end] - ff[1])
        Nii = @. (1.0 - fff) * Ni[1] + fff * Ni[end]
        hii = reverse(ip.(Nii))
        # << interpolate N, h

    end
    p = ip_hp.(h)
    T = ip_hT.(h)
    N = @. p / (c_kB * T)
    h_iout = ones(Int64, length(h))
    Atmosphere(h_iout, h, p, T, N)
end

#par = parameter()
#make_hpTN(par)