function PolyBairstov(a::Vector{Float64}, p::Float64=0.0, q::Float64=0.0, maxiter::Int64=1000, abstol::Float64=1e-8)
    n::Int64 = length(a)
    iters::Int64 = (isodd(n) ? (n - 3) / 2 : n / 2 - 1)
    g::Vector{Float64} = a[:]
    r::Float64 = Inf
    s::Float64 = Inf
    rq::Float64 = Inf
    sq::Float64 = Inf
    rp::Float64 = Inf
    sp::Float64 = Inf
    dp::Float64 = 0.0
    dq::Float64 = 0.0

    count::Int64 = 0

    res::Matrix{Float64} = Matrix{Float64}(undef, 2, Int64(iters + 1))

    for j = 1:iters
        count = 0
        g[1:n-2*j+2] = PolyDiv!(g[1:n-2*j+2], p, q)
        r = g[n-2*j+1]
        s = g[n-2*j+2]
        while (abs(r) + abs(s)) > abstol && count < maxiter
            g[1:n-2*j] = PolyDiv!(g[1:n-2*j], p, q)
            rq = -g[n-2*j-1]
            sq = -g[n-2*j]
            rp = sq - p * rq
            sp = -q * rq
            dq = (r * sp - s * rp) / (rp * sq - rq * sp)
            dp = (r * sq - s * rq) / (rp * sq - rq * sp)
            q += dq
            p -= dp
            count += 1

            g = a[:]
            g[1:n-2*j+2] = PolyDiv!(g[1:n-2*j+2], p, q)
            r = g[n-2*j+1]
            s = g[n-2*j+2]
        end
        a = g[:]
        res[1, j] = p
        res[2, j] = q
    end
    if isodd(n)
        res[1, iters+1] = g[2]/g[1]
        res[2, iters+1] = g[3]/g[1]
    else
        res[1, iters+1] = g[1]
        res[2, iters+1] = g[2]
    end
    return res
end

function PolyDiv!(g::Vector{T}, p::Float64, q::Float64) where {T<:Real}
    for i = 1:(length(g)-2)
        g[i+1] -= g[i] * p
        g[i+2] -= g[i] * q
    end
    return g
end