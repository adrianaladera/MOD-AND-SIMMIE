using LinearAlgebra

global const µ = 3.986e14  # m^3/s^2

# Helper function to take the column-wise norms of a matrix
norm_col(A) = [norm(A[:,i]) for i=1:size(A,2)]

function feval(x,p,u)
    # Function eval with input x formatted as:
    # r_1x, r_1y, r_1z, v_1x, v_1y, v_1z, r_2x, r_2y, ...
    # 3-vec of position followed by 3-vec of velocity for each node
    
    # Optionally modify the exponential decay constant (d) and maximum force amplitude (k) here
    d = 20000  # O(1e4 - 1e5) seems to be decent?
    k = 5  # Remember, gravity in LEO is ~ 9.8 m/s^2
    # Can also set these as arrays and choose k[i], d[i] in the controller func call in the nested for loop below
    
    N = length(x);
    f = zeros(eltype(x), N);
    
    Threads.@threads for i ∈ 1:6:N
        @views f[i:(i+2)] .= x[(i+3):(i+5)] # derivative of position is velocity
        positionsi = @view x[i:(i+2)]
        gravterm = -μ * norm(positionsi) ^ -3
        @views f[(i+3):(i+5)] .+= positionsi .* gravterm
        for j ∈ i:6:N
            if i != j
                positionsj = @view x[j:(j+2)]
                r_ij = positionsi - positionsj
                dist = norm(r_ij)
                accel = tanh_tol_controller(k,d,r_ij)
                @view(f[(i+3):i+5]) .+= accel
                @view(f[(j+3):j+5]) .-= accel
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

function tanh_tol_controller(k,d,r_ij)
    # Tanh controller, but below 1e-5 just returns 0
    dist = norm(r_ij)
    mag = (-k/2 * tanh(dist/d-5) + k/2)
    if mag < 1e-5
        return [0,0,0]
    end
    return mag * r_ij/dist
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
