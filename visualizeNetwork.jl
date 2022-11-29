using PlotlyJS

function visualizeNetwork(xs, fname="./example.html", skips=100)
    # skips : plot only every Nth time step to make the plot smaller
    N = floor(Int, length(xs[1,:])/6)
    
    trace = [PlotlyJS.scatter(x=xs[1:skips:end,1+j*6], y=xs[1:skips:end,2+j*6], z=xs[1:skips:end,3+j*6],mode="lines",type="scatter3d") for j in 0:N-1]
    p = PlotlyJS.plot(trace)
    open(fname, "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end
