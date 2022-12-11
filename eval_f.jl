using LinearAlgebra

global const µ = 3.986e14  # m^3/s^2

# Helper function to take the column-wise norms of a matrix
norm_col(A) = [norm(A[:,i]) for i=1:size(A,2)]

function feval(x,p,u)
    # Function eval with input x formatted as:
    # r_1x, r_1y, r_1z, v_1x, v_1y, v_1z, r_2x, r_2y, ...
    # 3-vec of position followed by 3-vec of velocity for each node
    
    # Optionally modify the exponential decay constant (d) and maximum force amplitude (k) here
    d = 25000  # O(1e4 - 1e5) seems to be decent?
    k = 5  # Remember, gravity in LEO is ~ 9.8 m/s^2
    # Can also set these as arrays and choose k[i], d[i] in the controller func call in the nested for loop below
    
    N = length(x);
    f = zeros(N);
    
    # Derivative of position is velocity
    f[1:6:end] = x[4:6:end]  # x axis
    f[2:6:end] = x[5:6:end]  # y axis
    f[3:6:end] = x[6:6:end]  # z axis
    
    # Calculate -µ/|r|^3 for each node
    # Reshape x into 3 x N/3 and take every other col -> get just positions, not velocities
    dists = norm_col(reshape(x, 3, floor(Int, N/3))[:,1:2:end])
    grav_terms = -µ * dists.^-3
    
    # Set gravitational acceleration as r * (-µ/|r|^3)
    f[4:6:end] = x[1:6:end].*grav_terms
    f[5:6:end] = x[2:6:end].*grav_terms
    f[6:6:end] = x[3:6:end].*grav_terms
    
    for i in 1:floor(Int, N/6)  # For each node,
        for j in 1:floor(Int, N/6)  # Loop over each other node
           if i != j  # And add the contributions from the interactions between nodes
                idx_i = (i-1)*6
                idx_j = (j-1)*6
                # Vector of position from node i to note j
                r_ij = x[1+idx_i:3+idx_i] - x[1+idx_j:3+idx_j]
                
                f[4+idx_i:6+idx_i] += tanh_controller(k,d,r_ij)
            end
        end
    end
    return f
end

function exp_controller(k,d,r_ij)
    # Original controller
    dist = norm(r_ij)
    return k*exp(-dist/d) * r_ij/dist
end

function sig_controller(k,d,r_ij)
    # Exponential controller
    dist = norm(r_ij)
    return k/(1+exp(dist/d-5)) * r_ij/dist
end

function tanh_controller(k,d,r_ij)
    # Tanh controller
    dist = norm(r_ij)
    return (-k/2 * tanh(dist/d-5) + k/2) * r_ij/dist
end

function energy_controller(k,d,r,v)
    # This one doesn't work as well maybe, but it conserves orbital energy
    # r, v are x[1+idx_i:3+idx_i], x[4+idx_i:6+idx_i]
    binormal = cross(r, v)
    return k/(1+exp(norm(r_ij)/d-5)) * binormal/norm(binormal)
end

function energy(x)
    # Calculate total specific kinetic + gravitational potential energy for each object
    # to check for energy conservation
    # Won't be entirely accurate because it doesn't account for potential between objects
    
    # gravitational: -mu/r
    # kinetic: 1/2 v^2
    
    dists = norm_col(reshape(x, 3, floor(Int, length(x)/3))[:,1:2:end])
    U = -µ ./dists
    
    vels = norm_col(reshape(x, 3, floor(Int, length(x)/3))[:,2:2:end])
    K = 1/2 * vels.^2
    
    return sum(U)+sum(K)
end
