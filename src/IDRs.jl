module IDRs

export idrs, syl, stein

# Copyright (c) 2015 M. Schauer, R. Astudillo,  M. B. van Gijzen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

function omega(t, s, angle = .7)
        ns = vecnorm(s)
        nt = vecnorm(t)
        ts = frodot(t,s)
        rho = abs(ts/(nt*ns))
        om = ts/(nt*nt)
        if rho < angle 
           om = om*angle/rho
        end
        return om
end
function frodot(A,B) 
         dot(A[:], B[:])
end

stein_op(X, A, B) = X + A*X*B
stein(A, B, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C)) = idrs(stein_op, (A,B), C,s,tol,maxit,X0)

syl_op(X, A, B) = A*X + X*B
syl(A, B, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C)) = idrs(syl_op, (A,B), C, s, tol, maxit, X0)


function idrs(op, args, C, s = 8, tol = 1E-8, maxit = length(C)^2, X0 = zeros(C))
    n = size(C,1)
    b = size(C,2)
    X = X0
    R = C - op(X, args...)
    tolc = tol*vecnorm(C)
    normR = vecnorm(R)
    RES = [normR]

    if normR <= tolc           # Initial guess is a good enough solution
       iter = 0                 
       return X0
    end

    P = rand(n,b,s)

    U = zeros(n,b,s)
    G = zeros(n,b,s)
    Q = zeros(n,b)
    V = zeros(n,b)
    M = eye(s,s)
    f = zeros(s,1)

    om = 1.
    iter = 0
    while normR > tolc && iter < maxit
        for i in 1:s,
            f[i] = frodot(R,P[:,:, i])
        end
        for k in 1:s 
        #
        # Solve small system and make v orthogonal to P:
            c = M[k:s,k:s]\f[k:s] 
            V = c[1]*G[:,:, k]
            
            Q = c[1]*U[:,:, k]
            for i = k+1:s
                for i1 in 1:n, j1 in 1:b
                    V[i1, j1] += c[i-k+1]*G[i1, j1, i]
                    Q[i1, j1] += c[i-k+1]*U[i1, j1, i]
                end
            end
        #
        # Compute new U[:,k] and G[:,k], G[:,k] is in space G_j
            V[:] = R - V
            U[:,:, k] = Q + om*V 
            G[:,:, k] = op(U[:,:, k], args...) 
        #
        # Bi-Orthogonalise the new basis vectors: 
           for i in 1:k-1
                alpha = frodot(P[:,:, i],G[:,:, k])/M[i,i]
                for i1 in 1:n, j1 in 1:b
                    G[i1, j1, k] -=  alpha*G[i1, j1, i]
                    U[i1, j1, k] -=  alpha*U[i1, j1, i]
                end

           end
        #
        # New column of M = P"*G  (first k-1 entries are zero)
            for i =k:s
                M[i,k] = frodot(G[:,:, k],P[:,:, i] )       
            end
        #
        #  Make r orthogonal to q_i, i = 1..k 
            beta = f[k]/M[k,k]
            R += - beta*G[:,:, k]
            X += + beta*U[:,:, k]
        
            normR = vecnorm(R)
            RES = [RES; normR ]
            iter = iter + 1
            if  normR < tolc || iter > maxit
                 return X
            end 
            if k < s 
                    f[k+1:s]   = f[k+1:s] - beta*M[k+1:s,k]
            end
        
        end
        #
        # Now we have sufficient vectors in G_j to compute residual in G_j+1
        # Note: r is already perpendicular to P so v = r
        #
        V[:] = R
        Q[:] = op(V, args...)
        om = omega(Q, R, 0.7)
        R[:] = R - om*Q
        X += om*V
        normR = vecnorm(R)
        iter = iter + 1
        RES = [RES; normR]
    end
    if normR > tol
        warn("idrs() did not converge")
    end
    return X
end

end # module