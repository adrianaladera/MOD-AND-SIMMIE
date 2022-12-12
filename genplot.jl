using LinearAlgebra, Distributions, Random, Rotations, ForwardDiff, SparseArrays
using Printf

global const µ = 3.986e14  # m^3/s^2

includet("eval_f.jl")
includet("forward_euler.jl")
includet("trap.jl")
includet("generate_random_initial_state.jl")
# includet("visualizeNetwork.jl")

f(x, p) = feval(x, p, 0)
f(x) = f(x, 0)

mean_a = 8e6
stddev_a = 1e2
mean_e = 0.01
stddev_e = 0.1

N = 20

x0 = generate_random_initial_state(N, mean_a, stddev_a, mean_e, stddev_e, true);

t0 = 0
ti = 1e4  # one orbit ~ 5e4 sec
num_steps = 100

# xs = euler(f, t0, ti, x0, num_steps);
xs = trapezoid(f, t0, ti, x0, 0, num_steps)

function plotorb(xs, n)
    mat = reduce(hcat, xs)
    p = plot(mat[1, :], mat[2, :], mat[3, :])
    for i ∈ 1:6:(6*n)
        plot!(p, mat[i, :], mat[i+1, :], mat[i+2, :])
    end
    return p
end
