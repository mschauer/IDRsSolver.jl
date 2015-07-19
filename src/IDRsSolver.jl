module IDRsSolver

export idrs, syl, stein

arand(C) = rand(size(C))
#adot(A,B) = vecdot(A, B)
adot(A,B) = dot(A[:], B[:])
anorm(A) = vecnorm(A)
axpy!(a, X, Y) = BLAS.axpy!(a, X, Y)

function omega(t, s, angle = .7)
    ns = anorm(s)
    nt = anorm(t)
    ts = adot(t,s)
    rho = abs(ts/(nt*ns))
    om = ts/(nt*nt)
    if rho < angle 
        om = om*angle/rho
    end
    om
end

linsys_op(X, A) = A*X
idrs(A, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C)) = idrs_core(linsys_op, (A), C,s,tol,maxit,X0)

stein_op(X, A, B) = X + A*X*B
stein(A, B, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C)) = idrs_core(stein_op, (A,B), C,s,tol,maxit,X0)

syl_op(X, A, B) = A*X + X*B
syl(A, B, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C)) = idrs_core(syl_op, (A,B), C, s, tol, maxit, X0)


function idrs_core{T}(op, args, C::T, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zero(C))

    X = X0
    R = C - op(X, args...)
    tolc = tol*anorm(C)
    normR = anorm(R)
    res = [normR]
	iter = 0                 

    if normR <= tolc           # Initial guess is a good enough solution
        return X0
    end
    
    Z = zero(C)

    P = T[arand(C) for k in 1:s]
    U = T[copy(Z) for k in 1:s]
    G = T[copy(Z) for k in 1:s]
    Q = copy(Z)
    V = copy(Z)
    
    M = eye(s,s)
    f = zeros(s)
    c = zeros(s)

    om = 1.
    iter = 0
    while normR > tolc && iter < maxit
        for i in 1:s,
            f[i] = adot(R,P[i])
        end
        for k in 1:s 

            # Solve small system and make v orthogonal to P

#            c = M[k:s,k:s]\f[k:s] 
            c = lufact!(M[k:s,k:s])\f[k:s]
            V[:] = G[k]
            scale!(c[1], V)
        
            Q[:] = U[k]
            scale!(c[1], Q)
            for i = k+1:s
                axpy!(c[i-k+1], G[i], V)
                axpy!(c[i-k+1], U[i], Q)
            end

            # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j

            #V = R - V
            scale!(-1., V)
            axpy!(1., R, V)
            
            U[k][:] = Q 
            axpy!(om, V, U[k])
            G[k] = op(U[k], args...) 

            # Bi-orthogonalise the new basis vectors

            for i in 1:k-1
                alpha = adot(P[i],G[k])/M[i,i]
                axpy!(-alpha, G[i], G[k])
                axpy!(-alpha, U[i], U[k])
            end

            # New column of M = P"*G  (first k-1 entries are zero)

            for i in k:s
                M[i,k] = adot(G[k],P[i])       
            end

            #  Make r orthogonal to q_i, i = 1..k 

            beta = f[k]/M[k,k]
            axpy!(-beta, G[k], R)
            axpy!(beta, U[k], X)
        
            normR = anorm(R)
            res = [res; normR]
            iter += 1
            if normR < tolc || iter > maxit
                return X
            end 
            if k < s 
                f[k+1:s] = f[k+1:s] - beta*M[k+1:s,k]
            end
        
        end

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r

        V[:] = R
        Q = op(V, args...)
        om::Float64 = omega(Q, R, 0.7)
#        R -= om*Q
#        X += om*V
        axpy!(-om, Q, R)
        axpy!(om, V, X)

        
        normR = anorm(R)
        iter += 1
        res = [res; normR]
    end
    if normR > tol
        warn("idrs() did not converge")
    end
    return X
end

end # module
