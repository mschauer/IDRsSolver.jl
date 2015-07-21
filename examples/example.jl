using IDRsSolver

println("FDM discretisation of a 3D convection-diffusion-reaction problem on a unit cube")
println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

# Define system

# Defaults:
h = 0.025
eps1 = 1.
beta = [0/sqrt(5), 250/sqrt(5), 500/sqrt(5)]
r = 400


# Generate matrix
m = round(Int, 1/h)-1
n = m*m*m
Sx = sparse(Tridiagonal((-eps1/h^2-beta[1]/(2*h))*ones(m-1),2*eps1/h^2*ones(m),(-eps1/h^2+beta[1]/(2*h))*ones(m-1)))
Sy = sparse(Tridiagonal((-eps1/h^2-beta[2]/(2*h))*ones(m-1),2*eps1/h^2*ones(m),-eps1/h^2+beta[2]/(2*h)*ones(m-1)))
Sz = sparse(Tridiagonal((-eps1/h^2-beta[3]/(2*h))*ones(m-1),2*eps1/h^2*ones(m),-eps1/h^2+beta[3]/(2*h)*ones(m-1)))
Is = speye(m)       
A = kron(kron(Is,Is),Sx) + kron(kron(Is,Sy),Is)+ kron(kron(Sz,Is),Is) -r*speye(n)

x = linspace(h,1-h,m)
sol = kron(kron(x.*(1-x),x.*(1-x)),x.*(1-x))
b = A*sol

println(" ");
println("The parameters of the problem are :")
println("Gridsize h = ", (h),";")
println("Number of equations = ", (n),";")
println("Diffusion parameter = ", (eps1),";")
println("Convection parameters = (",(beta[1]),",",(beta[2]),",",(beta[3]),");")
println("Reaction parameter = ",(r)," (Note: positive reaction parameter gives negative shift to matrix);")
println(" ")



# Defaults for the iterative solvers:

tol = 1e-8
maxiter = 1000
      
s = 4
println("IDR(4) iteration...")
time = @elapsed x, _ = idrs( A, b; s=s, tol=tol, maxiter=maxiter)
#       resvec = log10(resvec/resvec(1))

println("Final accuracy: ", norm(b-A*x)/norm(b))
println("CPU time: ", time, "s.")
      
