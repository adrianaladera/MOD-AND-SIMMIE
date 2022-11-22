using ForwardDiff, LinearAlgebra
function newton(f, x0, ϵf, ϵΔx, ϵxrel, maxiter = 256)
    xₖ = x0
    k = 1
    yₖ = f(xₖ)
    erryₖ = norm(yₖ, Inf)
    errΔxₖ = Inf
    relΔxₖ = Inf
    while k <= maxiter && (erryₖ > ϵf || errΔxₖ > ϵΔx || relΔxₖ > ϵxrel)
        println("iteration $k")
        Jf = ForwardDiff.jacobian(f, xₖ)[1]
        if length(Jf) == 1
            Jf = Jf[1]
        end
        x_next = xₖ .- Jf \ yₖ
        yₖ = f(x_next)
        erryₖ = norm(yₖ)
        errΔxₖ = norm(x_next .- xₖ)
        relΔxₖ = norm(x_next .- xₖ) / maximum(abs.(xₖ))
        xₖ = x_next
        k += 1
    end
    return xₖ
end

function newtonHCM(f, x0, ϵf, ϵΔx, ϵxrel; dq = 0.01)
    x = x0
    k = 1
    H(x, q) = (1-q) * (x - x0) + q * f(x)
    xs = Vector{Float64}(undef, length(0:dq:1))
    for q ∈ 0:dq:1
        x = newton(Base.Fix2(H, q), x, ϵf, ϵΔx, ϵxrel)
        xs[k] = x
        k += 1
        println("\nx($q) = $x\n")
    end
    return x, xs
end

function poissonf(ψ)
    [
        poissonhelper(-10., ψ[1], ψ[2]; N = length(ψ)),
        [poissonhelper(ψ[i - 1], ψ[i], ψ[i + 1]; N = length(ψ)) for i ∈ 2:(length(ψ) - 1)]...,
        poissonhelper(ψ[end - 1], ψ[end], 10.; N = length(ψ))
    ]
end

function poissonhelper(ψ₁, ψ₂, ψ₃; N)
    return -((-ψ₁ + 2 * ψ₂ - ψ₃) / ((1 / (N + 1))^2)) - (ℯ^(ψ₂) - ℯ^(-ψ₂))
end

function newtonkeepxs(f, x0, ϵf, ϵΔx, ϵxrel, maxiter = 256)
    xₖ = x0
    k = 1
    yₖ = f(xₖ)
    erryₖ = norm(yₖ, Inf)
    errΔxₖ = Inf
    relΔxₖ = Inf
    xs = x0
    while k <= maxiter && (erryₖ > ϵf || errΔxₖ > ϵΔx || relΔxₖ > ϵxrel)
        println("iteration $k")
        Jf = ForwardDiff.jacobian(f, xₖ)[1]
        if length(Jf) == 1
            Jf = Jf[1]
        end
        x_next = xₖ .- Jf \ yₖ
        xs = hcat(xs, x_next)
        yₖ = f(x_next)
        erryₖ = norm(yₖ)
        errΔxₖ = norm(x_next .- xₖ)
        relΔxₖ = norm(x_next .- xₖ) / maximum(abs.(xₖ))
        xₖ = x_next
        k += 1
    end
    return xₖ, xs
end

function geterrs(xₖ, xs)
    xs .-= xₖ
    return norm.(eachcol(xs), Inf)
end
