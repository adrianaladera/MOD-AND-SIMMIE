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

orbits = 5
t0 = 0
ti = 8e3 * orbits
num_steps = 5000 * orbits

xs = euler(f, t0, ti, x0, num_steps);

# before you run this, check that this isn't going to overwrite another file
using DelimitedFiles
writedlm("simulation_$(N)objects_$(orbits)orbits.txt", xs)
