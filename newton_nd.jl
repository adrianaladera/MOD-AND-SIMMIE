using LinearAlgebra
include("jac.jl")

function NewtonNd(evaluate_f, x0,p,u,errf,errDeltax,relDeltax,max_iter)
    k = 1;   
    X = Float64[];
    X = copy(x0);      
    f = evaluate_f(X[:,k], p, u);
    errf_k      = norm(f,Inf);
    errDeltax_k = Inf;
    relDeltax_k = Inf;

    while k<=max_iter && (errf_k>errf || errDeltax_k>errDeltax || relDeltax_k>relDeltax)
       Jf = Jac(evaluate_f,X[:,k], p, u);
       Deltax      = Jf\(-f); 
       X = cat(X, X[:,k] + Deltax, dims=(2,2)) 
       k = k+1;
       f = evaluate_f(X[:,k], p, u);
       errf_k = norm(f,Inf);
       errDeltax_k = norm(Deltax,Inf);
       relDeltax_k = norm(Deltax,Inf)/maximum(abs.(X[:,k]));
    end

    x = X[:,k];
    
    iterations = k-1; 
    if errf_k<=errf && errDeltax_k<=errDeltax && relDeltax_k<=relDeltax
       converged = 1;
#        if visualize
          @printf("Newton converged in %d iterations\n", iterations);
#        end
    else
       converged = 0;
       @printf("Newton did NOT converge! Maximum Number of Iterations reached\n");
    end
    
    
    return x
end
