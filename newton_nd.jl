using LinearAlgebra
include("jac.jl")

function NewtonNd(evaluate_f, x0,p,u,errf,errDeltax,relDeltax,max_iter)
    k = 1
    X = copy(x0)    
    f = evaluate_f(X[:,k], p, u)
    errf_k      = norm(f,Inf)
    errDeltax_k = Inf
    relDeltax_k = Inf

    while k<=max_iter && (errf_k>errf || errDeltax_k>errDeltax || relDeltax_k>relDeltax)
       Jf = Jac(evaluate_f,X[:,k], p, u)
       Deltax      = Jf\(-f) 
       X = cat(X, X[:,k] + Deltax, dims=(2,2)) 
       k = k+1
       f = evaluate_f(X[:,k], p, u)
       errf_k = norm(f,Inf)
       errDeltax_k = norm(Deltax,Inf)
       relDeltax_k = norm(Deltax,Inf)/maximum(abs.(X[:,k]))
    end

    x = X[:,k]
    
    iterations = k-1 
    if errf_k<=errf && errDeltax_k<=errDeltax && relDeltax_k<=relDeltax
       converged = 1
#        if visualize
          @printf("Newton converged in %d iterations\n", iterations)
#        end
    else
       converged = 0
       @printf("Newton did NOT converge! Maximum Number of Iterations reached\n")
    end
    return x
end

function newtonHCM(evaluate_f, x0,p,u,errf,errDeltax,relDeltax,max_iter; dq = 0.05)
   x = x0
   k = 1
   H(x, p, u, q) = (1-q) * (x - x0) + q * evaluate_f(x, p, u)
   for q ∈ 0:dq:1
      H2(x, p, u) = H(x, p, u, q)
      x = NewtonNd(H2, x, p, u, errf, errDeltax, relDeltax, max_iter)
      k += 1
      println("\nx($q) = $x\n")
   end
   return x
end