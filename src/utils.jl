using Interpolations


"""
    Interpolate data onto a equidistant grid with n grid points

    x0 Vector{Float64} -- [description]
    y0 Vector{Float64} -- [description]
    n Int64 -- [description]

    (x, y) -- [description]
"""
function lininterp(x0, y0, np)
    lip = linear_interpolation(x0, y0, extrapolation_bc = Line())
    x = collect(LinRange(x0[1], x0[end], np))
    y = lip.(x)
    x, y, lip
end
