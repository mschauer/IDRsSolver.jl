using IDRsSolver
using Base.Test

srand(51)
n = 10
maxiter = n*n*n
tol = 1E-8
s = 8

## real

A = -eye(n) + 0.7*(2rand(n,n)-1)
B = -eye(n) + 0.7*(2rand(n,n)-1)
xsol = rand(n)
Xsol = 0.5eye(n) + 0.3(2rand(n,n)-1)
C = Xsol + A*Xsol*B
D = A*Xsol + Xsol*B
b = A*xsol


X, _ = stein(A, B, C; s=s, tol=tol, maxiter=maxiter)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(C)

X, _ = syl(A, B, D; s=s, tol=tol, maxiter=maxiter)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(D)

x, _ = idrs(A, b; s=s, tol=tol, maxiter=maxiter)
@test norm(xsol - x) <= 10*tol*norm(b)

# complex

A = -eye(n) + 0.7*(2rand(n,n)-1) + im*0.7*(2rand(n,n)-1)
B = -eye(n) + 0.7*(2rand(n,n)-1) + im*0.7*(2rand(n,n)-1)
xsol = rand(n) + im*rand(n)
Xsol = 0.5eye(n) + 0.3*(2rand(n,n)-1) + im*0.3*(2rand(n,n)-1)
C = Xsol + A*Xsol*B
D = A*Xsol + Xsol*B
b = A*xsol


X, _ = stein(A, B, C; s=s, tol=tol, maxiter=maxiter)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(C)

X, _ = syl(A, B, D; s=s, tol=tol, maxiter=maxiter)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(D)

x, _ = idrs(A, b; s=s, tol=tol, maxiter=maxiter)
@test norm(xsol - x) <= 10*tol*norm(b)
