%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  1D Finite Element Method Solver 
%  This code solves the model problem -u''(x) = f for Omega = [x1,x2] with
%  Dirichlet boundary conditions u(x1) = u(x2) = 0
%   
%  This program solves the method for different values of h (step size) and
%  provides error analysis by producing the error for the infinity norm, 
%  L2 norm and H1 norm. 
%
%  x1, x2 - interval bounds
%  
%  node_no - number of nodes
%  element_no - number of elements
%  x - our array which contains the coordinates of our nodes
%  h - length of each element
%
%  [xi, w] - Our gaussian points and gaussian weights (respectively)
%
%  nodes - this will be our 2 x (element_no) array which contains the
%  numbering for each pair of nodes
%
%  A_global, B_global - This will be our global matrices for our stiffness
%  matrix and load vector, respectively
%
%  u - Our finite element solution
%
%
%  REFERENCE:  Numerical Solution of Differential Equations -- Introduction to 
%              Finite Difference and Finite Element Methods, Cambridge University Press, 2017



%%%%%%%%%%%%%%%%%%%%%%%%%      Main Section      %%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all

L2_errors = [];
H1_errors = [];
h_values = [];

fprintf('numElem\t\th\t\t\tL2Error\t\t\t\tH1Error\n');

for numNodes = [6 11 21 41 81 161 321 641]
    numElem = numNodes - 1;
    x1 = 0;
    x2 = 3;
    h = (x2 - x1)/numElem;
    [L2err, H1err, A_global] = main(numNodes);
    
    % Store errors and corresponding step sizes
    L2_errors = [L2_errors, L2err];
    H1_errors = [H1_errors, H1err];
    h_values = [h_values, h];
    fprintf('%d\t\t%.6f\t\t%.6e\t\t%.6e\n', numElem ,h, L2err, H1err);
end

% Plot errors vs. step sizes
loglog(h_values, L2_errors, 'ro-', h_values, H1_errors, 'bo-', 'LineWidth', 1.5)
hold on
loglog(h_values, 0.5*(L2_errors(1)/(h_values(1)^2))*h_values.^2, 'r--', ...
       h_values, 0.5*(H1_errors(1)/h_values(1))*h_values, 'b--', 'Linewidth', 1.5)
lgd = legend("L_2 error", "H_1 error", 'O(h^2)', 'O(h^1)', 'Location', 'southeast');
xlabel('Mesh size')  % x-axis label
ylabel('Error')  % y-axis label
title('Error Analysis for L2 and H1 norm')
lgd.FontSize = 10;


function [L2err, H1err, A_global] = main(node_no)

    % Define bounds of the domain
    x1 = 0;
    x2 = 3;
    
    % Define step size h
    element_no = node_no - 1;
    h = (x2-x1)/element_no;
    
    x = zeros(node_no,1);
    [xi, w] = quadrature();
    
    % Setup x vector for coordinates of nodes on the domain
    for i = 1:node_no
        x(i) = x1 + (i-1)*h;
    end
    
    % Setup matrix for storing the node pairs for each element
    nodes = zeros(2, element_no);
    for i = 1:element_no
        for j = 1:2
            nodes(j,i) = j + (i-1);
        end
    end
    
    % Initialise global stiffness matrix and load vector
    A_global = zeros(node_no, node_no);
    B_global = zeros(node_no, 1);
    
    % Assemble stiffness matrix and load vector
    [A_global, B_global] = form_matrices(nodes, x, element_no, xi, w, A_global, B_global, node_no);
    
    % Solve for u
    u = A_global\B_global;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%     Exact Solutions to u     %%%%%%%%%%%%%%%%%%%%%%%
    

    uexact = @(x) -(x-3)^2*x^2;
    uderivative = @(x) -4*x^3 + 18*x^2 - 18*x;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%      Calculate Errors     %%%%%%%%%%%%%%%%%%%%%%%%%
    
    [L2err, H1err] = error(uexact, uderivative, x, xi, w, nodes, u, element_no);
    
    x3 = linspace(x1, x2, 1000);
    u_analytical = u_exact(x3);
    
    % Plotting
    figure;
    plot(x, u, 's-', "LineWidth", 1.5);
    hold on;
    plot(x3, u_analytical, 'r-', 'LineWidth', 1.5);
    hold on;
    xlabel('x'); 
    ylabel('u(x)  vs.  u_{fem}(x)'); 
    title(sprintf('Finite Element Solution vs Exact Solution for numElem = %d', node_no-1));
    legend('Finite Element Solution', 'Exact Solution');
    hold off;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%      Post-Processing      %%%%%%%%%%%%%%%%%%%%%%%%%
    
    x3 = linspace(x1, x2, 1000);
    for i = 1:length(x3)
        u_analytical(i) = uexact(x3(i));
    end

end


%%%%%%%%%%%%     Numerical Quadrature and Reference Element    %%%%%%%%%%%%


function [xi, w] = quadrature()
    xi = [-1/sqrt(3); 1/sqrt(3)];
    w = [1; 1];
end

function [N, dN] = basis(xi)
    N = [(1-xi)/2; (1+xi)/2];
    dN = [-1/2; 1/2];
end


%%%%%%%%%%%%%%%%%%%%     Assemble Matrices A, B     %%%%%%%%%%%%%%%%%%%%%%%


function [A_global, B_global] = form_matrices(nodes, x, element_no, ...
    xi, w, A_global, B_global, node_no)
    for nel = 1:element_no
        i1 = nodes(1, nel);
        i2 = nodes(2, nel);

        [A_local, B_local] = local_mat(x(i1), x(i2), xi, w);
        [A_global, B_global] = global_mat(A_local, B_local, nel, nodes, A_global, B_global);
    end
    
    
    % Apply Dirichlet boundary conditions
    A_global(1, :) = 0;  % Zero out the first row
    A_global(:, 1) = 0;  % Zero out the first column
    A_global(1, 1) = 1;  % Set diagonal element to 1
    B_global(1) = 0;  % Set the value at the first node
    
    A_global(node_no, :) = 0;  % Zero out the last row
    A_global(:, node_no) = 0;  % Zero out the last column
    A_global(node_no, node_no) = 1;  % Set diagonal element to 1
    B_global(node_no) = 0;  % Set the value at the last node
end



%%%%%%%%%%%%%%%%%%%%      Calculate Local Matrices      %%%%%%%%%%%%%%%%%%%


function [A_local, B_local] = local_mat(x1, x2, xi, w)
    dx = (x2-x1)/2;
    xk = 1;


    A_local = zeros(2,2);
    B_local = zeros(2,1);

    for l=1:2
        x = x1 + (x2-x1)/2 * (1+xi(l));
        xf = f(x);
        [N,dN] = basis(xi(l));
        for i =1:2 
            B_local(i) = B_local(i) + N(i)*xf*w(l)*dx;
            for j =1:2 
                A_local(i,j) = A_local(i,j) + (xk*dN(i)*dN(j)/(dx*dx))*w(l)*dx;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%     Assemble Global Matrices   %%%%%%%%%%%%%%%%%%%%%


function [A_global, B_global] = global_mat(A_local, B_local, nel, nodes, A_global, B_global)
    for i = 1:2
        iglobal = nodes(i, nel);
        B_global(iglobal) = B_global(iglobal) + B_local(i);
        for j = 1:2
            jglobal = nodes(j, nel);
            A_global(iglobal, jglobal) = A_global(iglobal, jglobal) + A_local(i,j);
        end
    end
end





%%%%%%%%%%%%%%%%%%       Calculate the Error Norms       %%%%%%%%%%%%%%%%%%

function [L2err, H1err] = error(uexact, uderivative, x, xi, w, nodes, u, element_no)
    L2err = 0;
    H1err = 0;
    for nel = 1:element_no
        i1 = nodes(1, nel);
        i2 = nodes(2, nel);
        
        h = x(i2) - x(i1);

        errl2 = 0;
        errh1 = 0;
        for ig = 1:length(xi)
            [N, dN] = basis(xi(ig));
            x_interp = N(ig) * x(i1) + (1 - N(ig)) * x(i2); % Interpolated x within the element
            dx_interp = (dN(ig) * x(i1) + (1 - dN(ig)) * x(i2))/2 ; % Interpolated derivative of x within the element
            u_interp = N(ig) * u(i1) + (1 - N(ig)) * u(i2); % Interpolated u within the element



            u_deriv = uderivative(dx_interp);
            u_ex = uexact(x_interp); % Exact solution at the interpolated point
            u_deriv_interp = dN(ig) * u(i1) + (1 - dN(ig)) * u(i2);
            
            errl2 = errl2 + (u_ex - u_interp)^2 * w(ig);
            errh1 = errh1 + ((u_ex - u_interp)^2 + ((u_deriv - u_deriv_interp)*h).^2) * w(ig);
        end

        % Scale error by element size and accumulate
        L2err = L2err + errl2 * h;
        H1err = H1err + errh1*h;
    end

    % Take square root to obtain L2, H1 error norm
    L2err = sqrt(L2err);
    H1err = sqrt(H1err);
end

%%%%%%%%%%%%%     Exact Solution for u, values for f     %%%%%%%%%%%%%%%%

function y = u_exact(x)
    
    y = -(x-3).^2.*x.^2;

end

function y = f(x)
    y = 12*x^2-36*x+18;
end
