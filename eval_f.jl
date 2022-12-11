using LinearAlgebra

global const µ = 3.986e14  # m^3/s^2

# Helper function to take the column-wise norms of a matrix
norm_col(A) = [norm(A[:,i]) for i=1:size(A,2)]

function feval(x, p, u)
    d = 0.0001
    k = 100 #100
    
    N = length(x);
    f = zeros(eltype(x), N);

    Threads.@threads for i ∈ 1:6:N
        @views f[i:(i+2)] .= x[(i+3):(i+5)] # derivative of position is velocity
        positionsi = @view x[i:(i+2)]
        gravterm = -μ * norm(positionsi) ^ -3
        @views f[(i+3):(i+5)] .= positionsi .* gravterm
        for j ∈ i:6:N
            if i != j
                positionsj = @view x[j:(j+2)]
                r_ij = positionsi - positionsj
                dist = norm(r_ij)
                @view(f[(i+3):i+5]) .+= k * exp(-dist / d) * r_ij / dist
                @view(f[(j+3):j+5]) .-= k * exp(-dist / d) * r_ij / dist
            end
        end
    end
    return f
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
