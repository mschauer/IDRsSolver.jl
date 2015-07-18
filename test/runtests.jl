using IDRsSolver
using Base.Test

srand(55)
n = 10
b = 10
A = rand(n,n)
B = rand(b,n)
Xsol = rand(n,b)
C = Xsol + A*Xsol*B
D = A*Xsol + Xsol*B

maxit = (n*b)^2
tol = 1E-8
s = 8

X = stein(A,B, C, s, tol, maxit)
@test vecnorm(Xsol - X) <= 2tol*vecnorm(C)

X = syl(A,B, D, s, tol, maxit)
@test vecnorm(D - A*X-X*B) <= tol*vecnorm(D)
