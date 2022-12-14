using LinearAlgebra, Distributions, Random, Rotations
using Printf
using Plots

global const µ = 3.986e14  # m^3/s^2

include("eval_f.jl")
include("forward_euler.jl")
include("trapezoidal.jl")
include("generate_random_initial_state.jl")
include("visualizeNetwork.jl")

function f(x)
    feval(x,0,0)
end

mean_a = 1e7
stddev_a = 1e5
mean_e = 0.01
stddev_e = 0.2

N = 200

x0 = generate_random_initial_state(N, mean_a, stddev_a, mean_e, stddev_e, true);

t0 = 0
ti = 8e3 * 5
num_steps = 5000 * 5

xs = euler(f, t0, ti, x0, num_steps);

#plot(xs[:,1], xs[:,2], xs[:,3])

using DelimitedFiles
writedlm("orbs_large.txt", xs)
