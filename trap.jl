using LinearAlgebra, ForwardDiff, SparseArrays

f̂(x, p::P) where P = 
    x - p[1] - p[2] / 2 * p[3](x, p[end]) - p[2] / 2 * p[3](p[1], p[end])
function trapezoid(f::F, t₀, tₙ, x₀, p, n) where F
    Δ = (tₙ - t₀) / n
    t = zeros(n + 1)
    x = Vector{typeof(x₀)}(undef, n + 1)
    Jf = zeros(length(x0), length(x0))
    x[begin] = x₀
    for i ∈ 1:n
        # gross hack here with the tuple. I don't think this will infer well, so lots of
        # allocations and such...
        x[i + 1] = newton(f̂, x[i], (x[i], Δ, f, p), 1e-4, 1e-4, 1e-4, 50; Jf)[1]
    end
    return x
end

function newton(f::F, x, p, errf, errΔₓ, relΔₓ, itermax; Jf = zeros(length(x), length(x))) where F
    k = 1
    fₖ = f(x, p)
    errfₖ = norm(fₖ, Inf)
    errΔxₖ = Inf
    relΔxₖ = Inf
    fₚ(x) = f(x, p)
    while k <= itermax && (errfₖ > errf || errΔxₖ > errΔₓ || relΔxₖ > relΔₓ)
        Jf = ForwardDiff.jacobian!(Jf, fₚ, x)
        # display(SparseMatrixCSC(Jf))
        Jsparse = SparseMatrixCSC(Jf)
        Jsparse = droptol!(Jsparse,1e-5)
        #println(nnz(Jsparse))
        Δₓ = Jsparse \ (-fₖ)
        x += Δₓ
        k += 1
        fₖ = fₚ(x)
        errfₖ = norm(fₖ, Inf)
        errΔxₖ = norm(Δₓ, Inf)
        relΔxₖ = norm(Δₓ, Inf) / maximum(abs.(x))
    end
    return x, k-1, errfₖ ≤ errf && errΔxₖ ≤ errΔₓ && relΔxₖ ≤ relΔₓ
end
