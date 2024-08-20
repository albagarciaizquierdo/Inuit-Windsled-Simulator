function x_eq = equilibrium_conditions(X_red_eq)

x_eq       = zeros(18,1);

x_K_eq      = X_red_eq(1,1);
z_K_eq       = X_red_eq(2,1);
theta_K_eq   = X_red_eq(3,1);

x_eq(7,1)  = x_K_eq;
x_eq(9,1)  = z_K_eq;
x_eq(11,1) = theta_K_eq;

end

