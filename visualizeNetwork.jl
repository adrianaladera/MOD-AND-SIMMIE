using PlotlyJS

function visualizeNetwork(xs, fname="./example.html")
    N = floor(Int, length(xs[1,:])/6)
    
    trace = [PlotlyJS.scatter(x=xs[:,1+j*6], y=xs[:,2+j*6], z=xs[:,3+j*6],mode="lines",type="scatter3d") for j in 0:N-1]
    p = PlotlyJS.plot(trace)
    open(fname, "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end
