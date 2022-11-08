using LinearAlgebra

include("newton_nd.jl")

function trapezoidal_homotopy(f, x0, xn, y0, n)
    # Use the trapezoidal method to propagate the orbits
    # x0: initial time
    # xn: final time
    # n: number of steps
    # y0: initial state

    h = (xn - x0)/n  # time step
    xs = zeros(floor(Int,n+1))
    ys = zeros(floor(Int,n+1),length(y0))
    xs[1] = x0
    ys[1,:] = y0
    
    # step: nonlinear solve (for x1) of x1 - x0 - dt/2 f(x1) - dt/2 f(x0)
    function step(x1, x0, h)
    	return x1 - x0 - h/2 * f(x1) - h/2 * f(x0)
    end
    
    for i in 1:floor(Int,n)
        xs[i + 1] = xs[i] + h
        # For reference, this is what forward euler does at this step:
        #ys[i + 1,:] = ys[i,:] + h * f(ys[i,:])

        # Instead:
        ys[i+1,:] = newtonHCM(step, ys[i,:], ys[i,:], h, 1e-9, 1e-8, 1e-8, 10000)
        # (Initial guess for ys[i+1] is ys[i])
    end
    
    return ys
end

function trapezoidal(f, x0, xn, y0, n)
    # Use the trapezoidal method to propagate the orbits
    # x0: initial time
    # xn: final time
    # n: number of steps
    # y0: initial state

    h = (xn - x0)/n  # time step
    xs = zeros(floor(Int,n+1))
    ys = zeros(floor(Int,n+1),length(y0))
    xs[1] = x0
    ys[1,:] = y0
    
    # step: nonlinear solve (for x1) of x1 - x0 - dt/2 f(x1) - dt/2 f(x0)
    function step(x1, x0, h)
    	return x1 - x0 - h/2 * f(x1) - h/2 * f(x0)
    end
    
    for i in 1:floor(Int,n)
        xs[i + 1] = xs[i] + h
        # For reference, this is what forward euler does at this step:
        #ys[i + 1,:] = ys[i,:] + h * f(ys[i,:])

        # Instead:
        ys[i+1,:] = NewtonNd(step, ys[i,:], ys[i,:], h, 1e-9, 1e-8, 1e-8, 10000)
        # (Initial guess for ys[i+1] is ys[i])
    end
    
    return ys
end
