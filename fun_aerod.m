function [F_a,M_a,alpha,beta,v_a] = fun_aerod(p,x,v_OK,R)
%% Description: This function computes the aerodynamic force and moment on the kite
% Inputs: parameters, state vector & kite velocity vector
% Outputs: aerodynamic force, moment in S_K, AoA, sidelip angle (radians)

%% Vectors in S_E
i_E = [1; 0; 0];
j_E = [0; 1; 0];
k_E = [0; 0; 1];

%% Aerodynamic velocity
v_a = v_OK -(-p.v_w*i_E); % in S_E
v_a = R.R_KE*v_a; % in S_K

%% Angle of attack and sideslip angle
alpha = atan(v_a(3)/v_a(1)); % (rad)
beta = atan(v_a(2)/norm(v_a)); % (rad)

%% Aerodynamic force
F_a = 0.5*p.rho*p.S*(norm(v_a))^2*[p.Cx0+p.Cxalpha*alpha; p.Cybeta*beta; p.Cz0 + p.Czalpha*alpha]; % in S_K

%% Aerodynamic moment
P = x(16)*p.b/(2*p.V_T);    % (-)
Q = x(17)*p.b/p.V_T;    % (-)
R = x(18)*p.b/(2*p.V_T);    % (-)
M_a = 0.5*p.rho*p.S*(norm(v_a))^2*p.c*[p.eps*(p.Clbeta*beta+p.Clp*P); p.Cm0 + p.Cmalpha*alpha + p.Cmq*Q; p.eps*(p.Cnbeta*beta+p.Cnr*R )];  % in S_K
end