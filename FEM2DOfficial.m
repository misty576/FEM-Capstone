%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  2D Finite Element Method Solver 
%  This code solves the model problem -\nabla u = f for a sqaure (or rectangular)
%  domain \Omega = [xmin,xmax] x [ymin,ymax] (Can be modified to your liking).
%
%  xmin, xmax, ymin, ymax - the bounds of our domain
%  Nx, Ny - Number of elements in the x,y direction resp.
%
%  nele - number of elements. (We multiply by 2 since each 'square' is made up of 2 triangles)
%  
%  p(1,1), p(1,2), ... , p(1,n) <-- x coordinates of the nodal points
%  p(2,1), p(2,2), ... , p(2,n) <-- y coordinates of the nodal points
%
%  t(1,1), t(1,2), ... , t(1, nele) <-- index of the first node of an element
%  t(2,1), t(2,2), ... , t(2, nele) <-- index of the second node of the element
%  t(3,1), t(3,2), ... , t(3, nele) <-- index of the thrid node of the element
%
%  e(1,1), e(1,2), ... , e(1, nbc) <-- index of the beginning node of a boundary edge
%  e(2,1), e(2,2), ... , e(2, nbc) <-- index of the end node of a boundary edge
%
% 
%  REFERENCE:  Numerical Solution of Differential Equations -- Introduction to 
%              Finite Difference and Finite Element Methods, Cambridge University Press, 2017
%

close all
clear all


%%%%%%%%%%%%%%%%      Plotting FEM Solution and Error     %%%%%%%%%%%%%%%%%
errors = [];
meshSizes = [];

%%%%  State the type of domain here. ('square' or 'circle')  %%%%%%%%
type = "square";

if type == "square"
    fprintf('numElem (on each axis)\t\tInfNormError\t\tRate of Convergence\n');
elseif type == "circle"
    fprintf('h\t\t\t\t\tInfNormError\t\tRate of Convergence\n');
end

i = 1;
for numElem = [4,8,16,32,64]

    [ufem, u_exact, totalerr, p, e, t, err] = main(numElem,numElem, type);
    if type == "square"
        h(i) = 1/numElem;
    elseif type == "circle"
        h(i) = 2/numElem;
    end
    
    % Store errors and corresponding step sizes
    errors = [errors, totalerr];
    meshSizes = [meshSizes, numElem]; % Update here
    if i == 1
        if type == "square"
            fprintf('%d\t\t\t\t\t\t\t%.6e\t\n', numElem, totalerr);
        elseif type == "circle"
             fprintf('%f\t\t\t%.6e\t\n', 2/numElem, totalerr);
        end

    else
        rate = log2(errors(i)/errors(i-1));
        if type == "square"
            fprintf('%d\t\t\t\t\t\t\t%.6e\t\t%d\t\n', numElem, totalerr, abs(rate));
        elseif type == "circle"
            fprintf('%f\t\t\t%.6e\t\t%d\t\n', 2/numElem, totalerr, abs(rate));
        end
    end
    figure;
    subplot(1, 2, 1);
    pdemesh(p,e,t, 'XYData', ufem, 'ZData', ufem, 'Mesh','on','ColorMap','jet');
    
    if type == "square"
        title(['FEM Solution for h = ', num2str(1/numElem)]);
    elseif type == "circle"
        title(['FEM Solution for h = ', num2str(2/numElem)]);
    end
    subplot(1, 2, 2);
    pdemesh(p,e,t,'XYData', err, 'ZData', err, 'Mesh','on','ColorMap','jet');
    title('FEM error u_{exact} - u_{fem}');
    i = i+1;
end

% Plot errors vs. step sizes
figure;
loglog(h, errors, 'ro-', 'LineWidth', 1.5);
hold on
loglog(h, 0.5*(errors(1)/(h(1)^2))*h.^2, 'b--', 'Linewidth', 1.5)
lgd = legend("error", 'O(h^2)', 'Location', 'southeast');
xlabel('Mesh size')  % x-axis label
ylabel('Error')  % y-axis label
title('Error Analysis for Infinity Norm')
lgd.FontSize = 10;



function [ufem, u_exact, totalerr, p, e, t, err] = main(Nx,Ny, type)


%%%%%%%%%%%%%%%%%%      Set up Uniform Triangulation      %%%%%%%%%%%%%%%%

% Define the domain [xmin, xmax] x [ymin, ymax]
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;


% Generate nodal points for p
x = linspace(xmin, xmax, Nx+1);
y = linspace(ymin, ymax, Ny+1);
[X, Y] = meshgrid(x, y);
p = [X(:)'; Y(:)']; % This will be our coordinates for each nodal point

% Generate elements and assign to t
nele = Nx*Ny*2; % Total number of elements 
t = zeros(4, nele);
count = 1;
for i = 1:Nx
    for j = 1:Ny
        n1 = (j-1)*(Nx+1) + i;
        n2 = n1 + 1;
        n3 = n1 + Nx + 1;
        n4 = n3 + 1;

        % Assign nodes for indices of triangles
        t(:,count) = [n1; n2; n3; 1]; % First triangle
        t(:,count+1) = [n2; n4; n3; 1]; % Second triangle

        count = count + 2;
    end
end


% Generate boundary edges and assign to e
nbc = 2*(Nx + Ny); % Total number of boundary edges
e = zeros(2, nbc);
count = 1;

% Bottom boundary
for i = 1:Nx
    e(:,count) = [i; i+1];
    count = count + 1;
end

% Right boundary
for j = 1:Ny
    e(:,count) = [(j-1)*(Nx+1) + Nx + 1; j*(Nx+1) + Nx + 1];
    count = count + 1;
end

% Top boundary
for i = Nx:-1:1
    e(:,count) = [Ny*(Nx+1) + i + 1; Ny*(Nx+1) + i];
    count = count + 1;
end

% Left boundary
for j = Ny:-1:1
    e(:,count) = [j*(Nx+1)+1; (j-1)*(Nx+1)];
    count = count + 1;
end



%%%%%%%%%%%%%%%%%%%%%     Assembling and Solving     %%%%%%%%%%%%%%%%%%%%%%  


if type == "circle"
    [p, e, t] = initmesh("circleg",'hmax',1/(Nx/2));
end

% Extract the number of elements and nodes
[~, nelem] = size(t);
[~, nnode] = size(p);

% Initialize arrays to store nodal indices for each element
nodes = zeros(3, nelem);

% Loop through each element to extract nodal indices
for i = 1:nelem
    nodes(1,i) = t(1,i);
    nodes(2,i) = t(2,i);
    nodes(3,i) = t(3,i);
end

% Initialize global stiffness matrix and force vector
A_global = zeros(nnode, nnode);
B_global = zeros(nnode, 1);


% Define out quadrature points and weights 

% 3 point quadrature
%xi_eta = [0, 1/2; 1/2, 0; 1/2, 1/2]; % Quadrature points
%weights = [1/6, 1/6, 1/6]; % Weights for each quadrature point

% 1 point quadrature
%xi_eta = [1/3,1/3];
%weights = 1/2;

% 4 point quadrature
xi_eta = [1/3,1/3; 2/15, 11/15; 2/15, 2/15; 11/15, 2/15];
weights = [-27/96, 25/96, 25/96, 25/96];

% Loop through each element to build stiffness matrix and load vector using quadrature
for nel = 1:nelem
    % Get nodal indices and coordinates for the current element
    i1 = nodes(1,nel);
    i2 = nodes(2,nel);
    i3 = nodes(3,nel);

    x1 = p(1,i1); y1 = p(2,i1);
    x2 = p(1,i2); y2 = p(2,i2);
    x3 = p(1,i3); y3 = p(2,i3);


    % Calculate the area of the triangle (for the Jacobian)
    area = abs(det([1,x1,y1;1,x2,y2;1,x3,y3]))/2;
    
    % Initialize element stiffness matrix and force vector
    A_local = zeros(3, 3);
    B_local = zeros(3, 1);
    
    % Loop over quadrature points
    for q = 1:length(weights)

        % Transform (ξ,η) to (x,y) using the shape functions
        xi = xi_eta(q, 1);
        eta = xi_eta(q, 2);
        N = [(1 - xi - eta), xi, eta]; % Linear shape functions
        
        % Calculate derivatives of basis functions
        dNdxi = [-1, 1, 0;
                -1, 0, 1];
    
        % Calculation of the Jacobian matrix
        Jac = [x2 - x1, y2 - y1;
             x3 - x1, y3 - y1];
        

        % Calculation of derivatives of N with respect to x and y using chain rule
        dNdx_dy = 1/2 * (Jac \ dNdxi);
        
        dNdx = dNdx_dy(1, :);
        dNdy = dNdx_dy(2, :);

        x = N(1)*x1 + N(2)*x2 + N(3)*x3;
        y = N(1)*y1 + N(2)*y2 + N(3)*y3;
        
        % Compute the Jacobian for the transformation
        J = area*2; % Since the area of the master element is 1/2
        
        % Integrate stiffness matrix and force/load vector
        for i = 1:3
            B_local(i) = B_local(i) + N(i)*f(x,y)*J*weights(q); 

            for j = 1:3
                A_local(i,j) = A_local(i,j) + (dNdx(i)*dNdx(j) + dNdy(i)*dNdy(j))*J/(2);
            end
        end
    end
    
    % Assembly into global stiffness matrix and load vector
    for i= 1:3
        iglobal = nodes(i,nel);
        B_global(iglobal) = B_global(iglobal) + B_local(i);
        for j=1:3
            jglobal = nodes(j,nel);
            A_global(iglobal,jglobal) = A_global(iglobal,jglobal) + A_local(i,j);
        end
    end
end


% Extract the number of prescribed boundary conditions
[~, npres] = size(e);
g = zeros(npres, 1);
% Apply Dirichlet boundary conditions
for i = 1:npres
    xb = p(1, e(1,i)); 
    yb = p(2, e(1,i));
    
    g(i) = uexact(xb,yb);
end


% Modify global stiffness matrix and force vector for Dirichlet boundary conditions
for i = 1:npres
    nod = e(1,i);
    for k = 1:nnode
        B_global(k) = B_global(k) - A_global(k, nod)*g(i);
        A_global(nod,k) = 0;
        A_global(k,nod) = 0;
    end
    A_global(nod, nod) = 1;
    B_global(nod) = g(i);
end


% Solve equation to obtain FEM solution
ufem = A_global\B_global;



%%%%%%%%%%%%%%%%%%%%%%     Calculating Error     %%%%%%%%%%%%%%%%%%%%%%%%% 


% Obtain exact solution for u
u_exact = zeros(nnode, 1);
for i = 1:nnode
    xt = p(1,i); yt = p(2,i);
    u_exact(i) = uexact(xt, yt);
end
    % Calculate infinity norm error
    err = ufem - u_exact;
    totalerr = abs(norm(err,Inf));

end

%%%%%%%%%%%%%     Exact Solution for u, values for f     %%%%%%%%%%%%%%%%

function yp = uexact(x,y)
 
    % 1.
    %yp = (1-x*x - y*y)/(4);

    % 2. 
    %yp = -cos(pi*x);

    % 3. 
    % yp = -sin(pi*x)*cos(2*pi*y);

    % 4.
    %yp = x*x + y*y;

    % 5.
    yp = sin(pi*x)*sin(pi*y);

    % 6. 
    %yp = 1/4*(x^2+y^4)*sin(pi*x)*cos(4*pi*y);
return

end


function f = f(x,y)
    
    % 1.
    %f = 1;	
    
    % 2. 
    %f = -pi*pi*cos(pi*x);
    
    % 3.
    % f = -5*pi^2*cos(2*pi*y)*sin(pi*x);
    
    % 4. 
    %f = -4;
    
    % 5.
     f = 2*pi^2*sin(pi*x)*sin(pi*y);
    
    % 6.
    %f = 17*pi^2*cos(4*pi*y)*sin(pi*x)*(x^2/4 + y^4/4) - 3*y^2*cos(4*pi*y)*sin(pi*x) - (cos(4*pi*y)*sin(pi*x))/2 - x*pi*cos(pi*x)*cos(4*pi*y) + 8*y^3*pi*sin(pi*x)*sin(4*pi*y);
    
    return
end

