using LinearAlgebra, Distributions, Random, Rotations

global const µ = 3.986e14  # m^3/s^2


function generate_random_initial_state(N, mean_a, stddev_a, mean_ecc, stddev_ecc)
    # Generate a state vector x0 corresponding to N nodes (objects).

    # Each object is randomly generated. The final distribution will have the specified
    # mean and standard deviation of semimajor axis "a" (~"orbit size")
    # and eccentricity "ecc" (~"orbit shape": how circular vs elliptical it is).
    # The other 4 orbital parameters are [0,2pi] periodic and uniformly randomly distributed.
    
    
    # First define the distributions to sample from.
    # Acceptable semimajor axis range: [6.7e6 m, inf] (earth radius + 400 km buffer, infinity)
    # Acceptable eccentricity range: [0, 0.9] (perfectly circular, nearly escape trajectory)
    # Ecc anom (eccentric anomaly) is the only [0,2pi] uniform param related to orbital shape.
    a_distribution = truncated(Normal(mean_a, stddev_a), 6.7e6, Inf)
    ecc_distribution = truncated(Normal(mean_ecc, stddev_ecc), 0, 0.9)
    ecc_anom_distribution = Uniform(0,2*pi)

    # Then sample from the distributions N times each.
    a = rand(a_distribution, N)
    ecc = rand(ecc_distribution, N)
    ecc_anom = rand(ecc_anom_distribution, N)

    # Calculate true anomaly.
    tru_anom_y = sqrt.(1 .+ ecc) .* sin.(ecc_anom/2)
    tru_anom_x = sqrt.(1 .- ecc) .* cos.(ecc_anom/2)
    tru_anom = 2 .* [atan(tru_anom_y[i], tru_anom_x[i]) for i in 1:N]

    # Position as defined in the orbital plane:
    dist = a .* (ecc - 1 .* cos.(ecc_anom))
    r_xi = dist .* cos.(tru_anom)
    r_yi = dist .* sin.(tru_anom)
    r_zi = zeros(N)

    # Velocity as defined in the orbital plane:
    velo = sqrt.(µ.*a) / dist
    v_xi = velo .* (-sin.(ecc_anom))
    v_yi = velo .* (sqrt.(1 .- ecc.^2) .* cos.(ecc_anom))
    v_zi = zeros(N)
    
    # The code above calculated the shape of the orbits but not the orientation.
    # To get the orientation, apply a unique random rotation matrix to each object.
    rand_rotations = [rand(RotMatrix{3}) for i in 1:N]
    r = vec(reduce(vcat,[[r_xi[i], r_yi[i], r_zi[i]]' * rand_rotations[i] for i in 1:N]))
    v = vec(reduce(vcat,[[v_xi[i], v_yi[i], v_zi[i]]' * rand_rotations[i] for i in 1:N]))

    # Assemble the state vector from positions and velocities.
    x0 = zeros(N*6)
    x0[1:6:end] .= r[1:3:end]
    x0[2:6:end] .= r[2:3:end]
    x0[3:6:end] .= r[3:3:end]
    x0[4:6:end] .= v[1:3:end]
    x0[5:6:end] .= v[2:3:end]
    x0[6:6:end] .= v[3:3:end]

    return x0
end
