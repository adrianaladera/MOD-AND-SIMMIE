using FiniteDifferences

function Jac(f::F, x₀, p, u; Δ=1e-5) where F
    m = length(f(x₀, p, u))
    n = length(x₀)
    o = similar(x₀, m, n)
    onehot = Vector{Float64}(undef, n)
    for i ∈ 1:m # iterating rows, slow
        for j ∈ 1:n
            fill!(onehot, 0) # slow
            onehot[j] = 1
            o[i,j] = (f(x₀ .+ Δ .* onehot, p, u)[i] - f(x₀ .- Δ .* onehot, p, u)[i]) / (2Δ)
        end
    end
    return o
end
