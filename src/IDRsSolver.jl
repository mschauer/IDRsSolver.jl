module IDRsSolver

#### IterativeSolver.jl reporting protocol
export idrs, syl, stein, idrs_core

type ConvergenceHistory{T, R}
    isconverged::Bool
    threshold::T
    mvps::Int
    residuals::R
end

function empty!(ch::ConvergenceHistory)
    ch.isconverged = false
    ch.mvps = 0
    empty!(ch.residuals)
    ch
end

function push!(ch::ConvergenceHistory, resnorm::Number)
    push!(ch.residuals, resnorm)
    ch
end

####

# define abstract versions of rand, dot, norm and axpy!

arand!(C) = rand!(C)
#adot(A,B) = vecdot(A, B)
adot(A,B) = dot(A[:], B[:])
anorm(A) = vecnorm(A)
aaxpy!(a, X, Y) = BLAS.axpy!(a, X, Y)

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

linsys_op(x, A) = A*x
idrs(A, b, x0 = zeros(b); s = 8, tol = sqrt(eps(anorm(b))), maxiter = length(b)^2) = idrs_core(linsys_op, (A,), b, x0; s = s, tol=tol, maxiter=maxiter)

stein_op(X, A, B) = X + A*X*B
stein(A, B, C, X0 = zeros(C); s = 8, tol = sqrt(eps(anorm(C))), maxiter = length(C)^2) = idrs_core(stein_op, (A,B), C, X0; s = s, tol=tol, maxiter=maxiter)

syl_op(X, A, B) = A*X + X*B
syl(A, B, C, X0 = zeros(C); s = 8, tol = sqrt(eps(anorm(C))), maxiter = length(C)^2) = idrs_core(syl_op, (A,B), C, X0; s = s, tol=tol, maxiter=maxiter)


function idrs_core{T}(op, args, C::T, X0 = zero(C); s = 8, tol = sqrt(eps(anorm(C))), maxiter = length(C)^2)

    X = X0
    R = C - op(X, args...)
    tolc = tol 
    normR = anorm(R)
    res = [normR]
	iter = 0                 

    if normR <= tolc           # Initial guess is a good enough solution
        return X0, ConvergenceHistory(0<= res[end] < tolc, tolc, length(resnorms), resnorms)
    end
    
    Z = zero(C)

    P = T[arand!(copy(C)) for k in 1:s]
    U = T[copy(Z) for k in 1:s]
    G = T[copy(Z) for k in 1:s]
    Q = copy(Z)
    V = copy(Z)
    
    M = eye(s,s)
    f = zeros(s)
    c = zeros(s)

    om = 1.
    iter = 0
    while normR > tolc && iter < maxiter
        for i in 1:s,
            f[i] = adot(R,P[i])
        end
        for k in 1:s 

            # Solve small system and make v orthogonal to P

            c = \(LowerTriangular(M[k:s,k:s]),f[k:s])
            copy!(V, G[k])
            scale!(c[1], V)
        
            copy!(Q, U[k])
            scale!(c[1], Q)
            for i = k+1:s
                aaxpy!(c[i-k+1], G[i], V)
                aaxpy!(c[i-k+1], U[i], Q)
            end

            # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j

            #V = R - V
            scale!(-1., V)
            aaxpy!(1., R, V)
            
            copy!(U[k], Q)
            aaxpy!(om, V, U[k])
            G[k] = op(U[k], args...) 

            # Bi-orthogonalise the new basis vectors

            for i in 1:k-1
                alpha = adot(P[i],G[k])/M[i,i]
                aaxpy!(-alpha, G[i], G[k])
                aaxpy!(-alpha, U[i], U[k])
            end

            # New column of M = P"*G  (first k-1 entries are zero)

            for i in k:s
                M[i,k] = adot(G[k],P[i])       
            end

            #  Make r orthogonal to q_i, i = 1..k 

            beta = f[k]/M[k,k]
            aaxpy!(-beta, G[k], R)
            aaxpy!(beta, U[k], X)
        
            normR = anorm(R)
            res = [res; normR]
            iter += 1
            if normR < tolc || iter > maxiter
                return X, ConvergenceHistory(0<res[end]<tolc, tolc, length(res), res)
            end 
            if k < s 
                f[k+1:s] = f[k+1:s] - beta*M[k+1:s,k]
            end
        
        end

        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r

        copy!(V, R)
        Q = op(V, args...)
        om::Float64 = omega(Q, R, 0.7)
        aaxpy!(-om, Q, R)
        aaxpy!(om, V, X)

        
        normR = anorm(R)
        iter += 1
        res = [res; normR]
    end
    if normR > tolc
        warn("idrs() did not converge")
    end
    return X, ConvergenceHistory(0<res[end]<tolc, tolc, length(res), res)
end

end # module
