using LinearAlgebra, ForwardDiff, FiniteDifferences

function transientreference(
    f, t₀, tₙ, Δt, solver=forwardeuler, x₀=0., p=(;); kwargs...
)
    n = round(Int64, (tₙ - t₀) / Δt)
    solver(f, t₀, tₙ, x₀, p, n)
end

# p(x[i], u(t[i]), Δ, f, p_inner)


function trapezoid_finitediff(f::F, t₀, tₙ, x₀, p, n) where F
    Δ = (tₙ - t₀) / n
    t = zeros(n + 1)
    x = Vector{typeof(x₀)}(undef, n + 1)
    t[begin] = t₀
    x[begin] = x₀
    for i ∈ 1:n
        t[i + 1] = t[i] + Δ
        # gross hack here with the tuple. I don't think this will infer well, so lots of
        # allocations and such...
        x[i + 1] = newton_finitediff(f̂, x[i], (x[i], Δ, f, p), t[i], 1e-4, 1e-4, 1e-4, 100)[1]
    end
    return x, t
end

"""
    forwardeuler(f(x, p, t), x₀, xₙ, y₀, p, t, n)
"""
function forwardeuler(f::F, t₀, tₙ, x₀, p, n) where F
    Δ = (tₙ - t₀) / n
    ts = zeros(n + 1)
    xs = Vector{typeof(x₀)}(undef, n + 1)
    ts[begin] = t₀
    xs[begin] = x₀
    for i ∈ 1:n
        ts[i + 1] = ts[i] + Δ
        xs[i + 1] = xs[i] + Δ * f(xs[i], p, ts[i])
    end
    return xs, ts
end

"""
    newton(f(x, p, t), x0, u, p, errf, errΔₓ, relΔₓ, itermax, Jf = nothing)

Use Newton's Method to solve the nonlinear system f(x, p, t) = 0

"""
function newton_finitediff(f::F, x, p, t, errf, errΔₓ, relΔₓ, itermax) where F
    k = 1
    X = [x]
    fₖ = f(X[k], p, u)
    errfₖ = norm(fₖ, Inf)
    errΔxₖ = Inf
    relΔxₖ = Inf
    fₚᵤ(x) = f(x, p, t)
    while k <= itermax && (errfₖ > errf || errΔxₖ > errΔₓ || relΔxₖ > relΔₓ)
        Jf = jacobian(central_fdm(5, 1), fₚᵤ, X[k])[1]
        if length(Jf) == 1
            Jf = Jf[1]
        end
        Δₓ = Jf \ (-fₖ)
        push!(X, X[k] + Δₓ)
        k += 1
        fₖ = fₚᵤ(X[k])
        errfₖ = norm(fₖ, Inf)
        errΔxₖ = norm(Δₓ, Inf)
        relΔxₖ = norm(Δₓ, Inf) / maximum(abs.(X[k]))
    end
    return X[k], k-1, X, errfₖ ≤ errf && errΔxₖ ≤ errΔₓ && relΔxₖ ≤ relΔₓ
end

H(f::F, x, p, t, q, x0) where F = (1-q) * (x - x0) + q * f(x, p, t)
"""
    newtonHCM(f(x, p, t), h(f, x, p, t, q, x0),  x0, t, p, ϵf, ϵΔₓ, ϵxrel; dq=0.01)
"""
function newtonHCM(f::F, h::H, x0, p, t, errf, errΔₓ, relΔₓ; dq = 0.01) where {F<:Function, H<:Function}
    X = [x0]
    k = 1
    q = Float64[0]
    converged = false
    while q[k] < 1
        x0 = X[k]
        q[k] + dq > 1 ? (push!(q, 1)) : push!(q, q[k] + dq)

        x, iterations, _, converged = newton((x, p, t)->h(f, x, p, t, q[k], x0), x0, p, t, errf, errΔₓ, relΔₓ, 1000)
        if !converged
            push!(q, q[k] - dq)
            dq = dq / 2
        else
            k = k + 1
            push!(X, x)
            if iterations <= 5
                dq = dq * 2
            else
                dq = dq * 1.25
            end
        end
    end
    return X[end], k, X, converged
end

