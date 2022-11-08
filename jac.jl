using FiniteDifferences

function Jac(f::F, x₀; Δ=1e-5) where F
    m = length(f(x₀))
    n = length(x₀)
    o = similar(x₀, m, n)
    for i ∈ 1:m
        for j ∈ 1:n
            onehot = zeros(n)
            onehot[j] = 1
            o[i,j] = (f(x₀ .+ Δ .* onehot)[i] - f(x₀ .- Δ .* onehot)[i]) / (2Δ)
        end
    end
    return o
end