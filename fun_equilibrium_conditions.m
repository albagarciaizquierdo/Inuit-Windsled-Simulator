function x_eq = fun_equilibrium_conditions(X_red_eq)
%% Description: 
% This function gives the full equilibrium state vector x_eq from the reduced
% equilibrium state vector X_red_eq
%% Inputs: 
% X_red_eq --> reduced equilibrium state vector (with x_K, z_K, theta)
%% Outputs: 
% x_eq --> equilibrium state vector

%%

x_eq       = zeros(18,1);

x_K_eq      = X_red_eq(1,1);
z_K_eq       = X_red_eq(2,1);
theta_K_eq   = X_red_eq(3,1);

x_eq(7,1)  = x_K_eq;
x_eq(9,1)  = z_K_eq;
x_eq(11,1) = theta_K_eq;

end

