function make_reduced_order(A, b, c, q)
    lambdas = eigvals(A);
    vectors = eigvecs(A);
    vectors
    Vq = zeros(N, q);
#     n = size(vectors)
#     print(n)
    for i in 0:q-1
        print(q-i)
        Vq[:,i+1] = vectors[:,q-i]
    end    
    A_hat =  Vq' * A * Vq;
    b_hat = Vq' * b;
    c_hat = Vq' * c;
#     a = size(Vq)
    return A_hat, b_hat, c_hat, Vq
end