# FEM-Capstone
The following is my repository for the code used in the examples for my Final Year Project titled "Finite Element Methods with Applications in Mathematical Finance"


DESCRIPTION OF FILES

1. FEM1DOFFICIAL.m 
- This MATLAB program solves the 1D Poisson Equation with Dirichlet Boundary Conditions and returns solutions for different sizes of h (= length of each element)
- The program incorporates Gaussian Quadrature (reference element, gaussian points and weights) to perform numerical integration.
- Performs error analysis by calculating the error of the FEM approximation wrt. the L2 and H1 error norm

2. FEM2DOFFICIAL.m 
- This MATLAB program solves the 2D Poisson Equation with Dirichlet Boundary Conditions and returns solutions for different levels of discretisation.
- The program incorporates Gaussian Quadrature (reference element, gaussian points and weights) to perform numerical integration.
- Performs error analysis by calculating the error of the FEM approximation wrt. the infinity error norm
- Can also be used to solve Poisson equations on non-conventional boundaries (circle or L-shaped domain)

3. BSFEMCall.m 
- This MATLAB program solves the Black-Scholes Equation (transformed to a heat equation) for the European Call option problem with given parameters.
- Performs error analysis by calculating and plotting the error of the FEM approximation

4. BSFEMPut.m 
- This MATLAB program solves the Black-Scholes Equation (transformed to a heat equation) for the European Put option problem with given parameters.
- Performs error analysis by calculating and plotting the error of the FEM approximation
