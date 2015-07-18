using IDRsSolver
using Base.Test

n = 10
b = 10
A = rand(n,n)
B = rand(b,n)
Xsol = rand(n,b)
C = Xsol + A*Xsol*B
X0 = rand(n,b)
maxit = n*n*n
tol = 1E-10
s = 10

X = stein(A,B, C,s,tol, maxit)
@test vecnorm(Xsol - X) <= tol
