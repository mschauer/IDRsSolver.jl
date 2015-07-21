IDRsSolver.jl
============================

The Induced Dimension Reduction method is a family of simple and fast Krylov
subspace algorithms for solving large nonsymmetric linear systems. The idea
behind the IDR(s) variant is to generate residuals that are in the nested
subspaces of shrinking dimension s. 

The function idrs() solves a general linear matrix equation

                         0 = C - op(X, args...) 
 
where op is a linear operator in X, for example 
        
        X -> X + A*X*B (Stein equation)

or

        X -> A*X + X*B (Sylvester equation).
        
                

Syntax
------

        X, h::ConvergenceHistory = idrs_core{T}(op, args, C::T, X0 = zero(C); s = 8, tol = sqrt(eps(anorm(C))), maxiter = length(C)^2)


        The right hand side C must support vecnorm, vecdot, copy!, rand! and axpy! or import and overload
        the abstract functions adot, anorm, arand! and aaxpy! from module IDRsSolver used by the procedure.
        .

Synonyms
--------
        idrs(A, b, ...) = idrs_core((x,A) -> A*x, (A,), b, ...) solves the linear equation equation Ax = b.
        stein(A, B, C, ...) = idrs_core((X,A,B) -> X + A*X*B, (A, B), C, ...) solves the Stein equation.
        syl(A, B, C, ...) = idrs_core((X,A,B) -> A*X + X*B, (A, B), C, ...) solves the Sylvester equation.

Arguments
---------
       
       s -- dimension reduction number. Normally, a higher s gives faster convergence, 
            but also  makes the method more expensive.
       tol -- tolerance of the method.  
       maxiter -- maximum number of iterations

       x0 -- Initial guess.
       
Output
------

        X -- Approximated solution by IDR(s)
        h -- Convergence history
        

        The [`ConvergenceHistory`](https://github.com/JuliaLang/IterativeSolvers.jl/issues/6) type provides information about the iteration history. 
            - `isconverged::Bool`, a flag for whether or not the algorithm is converged.
            - `threshold`, the convergence threshold
            - `residuals::Vector`, the value of the convergence criteria at each iteration        
        
        
        If ||C-op(X)||_F > tol, the function gives a warning.

    
References
----------

    [1] IDR(s): a family of simple and fast algorithms for solving large 
        nonsymmetric linear systems. P. Sonneveld and M. B. van Gijzen
        SIAM J. Sci. Comput. Vol. 31, No. 2, pp. 1035--1062, 2008 
    [2] Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits 
        Bi-orthogonality Properties. M. B. van Gijzen and P. Sonneveld
        ACM Trans. Math. Software,, Vol. 38, No. 1, pp. 5:1-5:19, 2011
    [3] This file is a translation of the following MATLAB implementation:
        http://ta.twi.tudelft.nl/nw/users/gijzen/idrs.m
    [4] IDR(s)' webpage http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html
