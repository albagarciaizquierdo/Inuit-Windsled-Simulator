function f  = fun_equilibrio_red(t,x_red)

x_K       = x_red(1,1);
z_K       = x_red(2,1);
theta_K   = x_red(3,1);

x_s       = zeros(18,1);
x_s(7,1)  = x_K;
x_s(9,1)  = z_K;
x_s(11,1) = theta_K;

dxdt = RHS(0,x_s);

f(1,1) = dxdt(13,1);
f(2,1) = dxdt(15,1);
f(3,1) = dxdt(17,1);

end