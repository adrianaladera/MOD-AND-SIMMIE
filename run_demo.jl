using DelimitedFiles

global const µ = 3.986e14  # m^3/s^2

include("eval_f.jl")
include("forward_euler.jl")
include("generate_random_initial_state.jl")
include("visualizeNetwork.jl")

f(x) = feval(x,0,0)

function run_sim(N, orbits)
    mean_a = 8.5e6
    stddev_a = 1e5
    mean_e = 0.01
    stddev_e = 0.1

    global x0 = generate_random_initial_state(N, mean_a, stddev_a, mean_e, stddev_e, true)

    for i in 1:orbits
        t0 = 0
        ti = 8e3
        num_steps = 5000
        xs = euler(f, t0, ti, x0, num_steps)
        global x0 = xs[end,:]
        writedlm("simulation_demo_$(i).txt", Float32.(xs[1:40:end-1,:]))
        println("completed orbit $(i)")
    end
    xs = vcat([readdlm("simulation_demo_$(i).txt") for i in 1:orbits]...)
    visualizeNetwork(xs,"./example.html",10)
    println("")
    println("saved visualization")
end
