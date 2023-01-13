#using DataInterpolations
#using DataFrames

#mutable struct LinInterpolator
#    linint :: Vector{LinearInterpolation}
#end
#
#function y_at(ip::LinInterpolator, x::Float64; ic=1)
#    ip.linint[ic](x)
#end
#
#function LinInterpolator(fname; header=false)
#    xy = read_npy(fname)
#    n1, n2 = size(xy)
#    x = xy[:,1]
#    intpols = Vector{LinearInterpolation}(undef,0)
#    for i in 2:n2
#        y = xy[:,i]
#        push!(intpols, LinearInterpolation(y,x))
#    end
#    LinInterpolator(intpols)
#end

mutable struct Interpolator
    xy  :: Matrix{Float64}
    x0  :: Float64
    ddx :: Float64
    nc  :: Int64
end

function Interpolator(fname)
    xy = read_npy(fname)
    nr, nc = size(xy)
    x0 = xy[1,1]
    ddx = 1.0 / (xy[1,2] - xy[1,1])
    Interpolator(xy, x0, ddx, nc)
end

function y_at(ip::Interpolator, x::Float64; ir=1)
    ix = max(1, min(floor(Int64, ((x - ip.x0) * ip.ddx))+1, ip.nc-1))
    u = (x - ip.xy[1,ix]) * ip.ddx
    (1.0 - u) * ip.xy[ir+1,ix] + u * ip.xy[ir+1,ix+1]
end
