%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 1D Black-Scholes Equation Solver for
%                        European Call Option
%
%  The following code solves the B-S equation as a transformed heat
%  equation using the Finite Element Method. This program solves the 
%  equation for different values of N and M (no. of spatial and temporal 
%  points resp.) and plots the solution. We also perform basic error 
%  analysis by taking the analytical solution and solving for the error.
%
%  T - Time to Maturity
%  K - Strike Price
%  r - risk free interest rate
%  sigma - volatility
%
%  x_min, x_max - stock price bounds (these are implemented so we can
%  define the spatial domain)
%  x, tau - space and time domain resp.
%  x_eff - spatial domain excluding boundaries
%
%  S, k1, r - variables for transformation to heat equation
%
%  b - matric used to describe the boundary conditions as described in
%  Section 5
%
%  xi - matrix for the time variables (calculated using Crank-Nicolson)
%
%  REFERENCE: The Finite Element Method for Option Pricing under Heston’s Model.
%             Final Report for MA6621 Programming and Computing for Finance and Actuarial Science.
%             Hu Wei, Yao Wuguannan and Huang Jiaheng, Univeristy of Hong Kong, 2016
%
close all
clear all


%%%%%%%%%%%%%%%%%%%%   Initial Setup of Variables   %%%%%%%%%%%%%%%%%%%%%%%

T = 0.5; K = 120;
r = 0.02; sigma = 0.35;


% Log transformation of stock price bounds
x_min = log(20/K); x_max = log(1500/K);



% Discretise the spatial and time domain

N = 100; M = 50;

x = linspace(x_min, x_max, N+2); dx = x(2) - x(1);
tau = linspace(0, T*(sigma^2)/2, M+1); dtau = tau(2) - tau(1);
x_eff = x(2 : length(x)-1);

xi = zeros(length(x_eff), length(tau));
b=xi;

k1 = 2*r / (sigma^2);
S = K*exp(x);


%%%%%%%%%%%%%%%%%%%  Establish boundary conditions  %%%%%%%%%%%%%%%%%%%%%%%


% Initial condition for option value at time T (this is gamma in our notes)
IC = @(x, q) max(0, exp(0.5 * x * (q + 1)) - exp(0.5 * x * (q - 1)));

% Boundary conditions (alpha - value at x_min , beta - value at x_max)
beta = @(tau) exp(1/2*(k1+1)*x_max + 1/4*tau*(k1+1)^2);
dbeta = @(tau) (1/4*tau*(k1+1)^2)*beta(tau);
alpha = @(tau) 0;
phi_boundary = @(x, tau) (beta(tau) - alpha(tau)) .* (x - x_min)/(x_max - x_min) + alpha(tau);

% Set initial condition.

xi(:,1) = IC(x_eff, k1)-phi_boundary(x_eff,0);

% Construct b to incorporate the boundary conditions into the FEM solution.
for i = 1 : length(tau)
    b(:, i) = ((x_eff - x_min)/(x_max - x_min)) * dx * (0.25 * (k1 + 1) ^ 2) * beta(tau(i));
end


%%%%%%  Assemble matrices and solve for xi using time stepping scheme  %%%%


A1 = zeros(N,N);
for i=1:N
    if i > 1
        A1(i, i-1) = -1/dx;
    end
    A1(i, i) = 2/dx;
    if i < N
        A1(i, i+1) = -1/dx;
    end
end
A1(1,1:2)=[2/dx,-1/dx]; A1(N,N-1:N)=[-1/dx,2/dx];

B1 = zeros(N,N);
for i=1:N
    if i > 1
        B1(i, i-1) = dx/6;
    end
    B1(i, i) = 2*dx/3;
    if i < N
        B1(i, i+1) = dx/6;
    end
end
B1(1,1:2)=[2*dx/3,dx/6]; B1(N,N-1:N)=[dx/6,2*dx/3];


% Solve for xi using the Crank-Nicolson Scheme
B_final = B1 + (0.5 * dtau) .* A1;
A_final = B1 - (0.5 * dtau) .* A1;

for i = 2 : length(tau)
    xi(:, i) = B_final \ (A_final * xi(:, i - 1) - (dtau/2) * (b(:, i) + b(:, i-1)));
end

%%%%%%%%%%%%%%%%%%%   Option Price Calculation V(S,t)   %%%%%%%%%%%%%%%%%%%

% Get option price at nodes
base = zeros(length(x), length(tau));
for i = 1 : length(x)
    for j = 1 : length(tau)
        base(i, j) = phi_boundary(x(i), tau(j));
    end
end

nodes = [zeros(1, length(tau)); xi; zeros(1, length(tau))];
nodes = nodes + base;


% Transform option prices back to the stock price and time domain
for i=1:length(x)
    for j=1:length(tau)
        nodes(i,j)=(K*exp((-0.5)*(k1-1)*x(i)+(-0.25)*((k1+1)^2)*tau(j)))*nodes(i,j);
    end
end



%%%%%%%%%%%%%%%%%%%%%%      Post-Processing      %%%%%%%%%%%%%%%%%%%%%%%%%%



% Output a 3D mesh of our FEM approiximation
V = nodes;
stock = S;
time = T-tau*2/sigma^2;
s = mesh(time',stock(N/8:N/2),V(N/8:N/2,:));view(3);
s.FaceColor = 'interp';
s.FaceAlpha = 0.3;
s.EdgeAlpha = 1; 
s.LineWidth = 0.75;
colormap turbo

%mesh(time', stock, Nodes_v);view(3);
title('Black-Scholes FEM Approximation for European Call');
ylabel('Stock Price (S)');xlabel('Time (T-t)');zlabel('Option Price V(S,t)');


% Output a 2D plot of the FEM Approximation for different time values
figure(2);
plot(stock(N/8:N/2), V(N/8:N/2, 1), 'r-', 'LineWidth', 1);
hold on;
plot(stock(N/8:N/2), V(N/8:N/2, M/2), '-', 'LineWidth', 1, 'Color',[0.1,0.6,0]);
plot(stock(N/8:N/2), V(N/8:N/2,M), 'b-', 'LineWidth', 1);
legend("V(S,T) at maturity", "V(S,T/2)", "V(S,0) initial condition");
xlabel('Stock Price (S)');
ylabel('Option Value V(S,t)');
title('European Call Option Prices at Different Times');



%%%%%%%%%%%%%%%%%%%%%%%      Error Analysis      %%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the analytical solution to the B-S Equation (using financial toolbox for MATLAB)
for i = 1:length(x)
    for j = 1:length(tau)
        [Call, Put] = blsprice(K*exp(x(i)), K, r, T - (tau(j)*2)/sigma^2, sigma);
        V_exact(i,j) = Call;
    end
end



V_exact = flip(V_exact,2);
V_exact = V_exact(1:3*N/4, :);
V = V(1:3*N/4,:);
err = V - V_exact;
figure(3);
mesh(time, stock(1:3*N/4), err);
colormap hsv
title('Error of Finite Element Solution');
ylabel('Stock Price');xlabel('Time');zlabel('Error');


