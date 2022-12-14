using LinearAlgebra

function euler(f, x0, xn, y0, n)
    # Use the forward euler method to propagate the orbits
    # x0: initial time
    # xn: final time
    # n: number of steps
    # y0: initial state

    h = (xn - x0)/n  # time step
    xs = zeros(floor(Int,n+1))
    ys = zeros(floor(Int,n+1),length(y0))
    xs[1] = x0
    ys[1,:] = y0
    
    for i in 1:floor(Int,n)
        xs[i + 1] = xs[i] + h
        ys[i + 1,:] = ys[i,:] + h * f(ys[i,:])
        if i%2000==0
            println(i)
        end
    end
    
    return ys
end
