using IDRsSolver
using Base.Test

srand(51)
n = 10
A = -eye(n) + 0.7*(2rand(n,n)-1)
B = -eye(n) + 0.7*(2rand(n,n)-1)
xsol = rand(n)
Xsol = 0.5eye(n) + 0.3(2rand(n,n)-1)
C = Xsol + A*Xsol*B
D = A*Xsol + Xsol*B
b = A*xsol

maxit = n*n*n
tol = 1E-8
s = 8

X = stein(A, B, C, s, tol, maxit)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(C)

X = syl(A, B, D, s, tol, maxit)
@test vecnorm(Xsol - X) <= 10*tol*vecnorm(D)

x = idrs(A, b, s, tol, maxit)
@test norm(xsol - x) <= 10*tol*norm(b)
