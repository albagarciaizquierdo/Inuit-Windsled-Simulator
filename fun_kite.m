function [F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a] = fun_kite(p,x,R,OK_AKp,OK_AKm,T_ASp,T_ASm)
%% Description: 
% This function computes forces and moments acting ond the
% kite
%% Inputs:
% x --> state vector/matrix
% p --> struct containing parameters
% R --> struct containing rotation matrices
% OK_AKp --> Position vector from OK to AK+
% OK_AKm --> Position vector from OK to AK-
% T_ASp --> Tether tension from AS+
% T_ASm --> Tether tension from AS-

%% Outputs: (all in S_K in SI units)
% F_K --> Kite force
% M_OK --> Kite moment about OK
% v_OK --> Kite velocity
% omega_KE --> Kite angular velocity
% H_OK --> Kite angular momentum
% W_K --> Kite weight
% T_AKp --> Tether tension from AK+
% T_AKm --> Tether tension from AK-
% F_a --> Aerodynamic force
% M_a --> Aerodynamic moment
% alpha --> angle of attack
% beta --> sideslip angle
% v_a --> Aerodynamic velocity

%% Vectors in S_E
i_E = [1; 0; 0];
j_E = [0; 1; 0];
k_E = [0; 0; 1];

%% Velocities and intertia
% Kite velocity
v_OK = [x(13); x(14); x(15)]; % in S_K
% v_OK = R_EK*v_OK; % in S_E
% Angular velocity omega_KE
omega_KE = [x(16); x(17); x(18)]; % in S_K;
% Inertia
I_OK = [p.Ix_K 0 0; 0 p.Iy_K 0; 0 0 p.Iz_K]; % in S_K
% Angular momentum 
H_OK = I_OK*omega_KE; % in S_K

%% Forces & moments
% Aerodynamic force and moment
[F_a,M_a,alpha,beta,v_a] = fun_aerod(p,x,v_OK,R); % in S_K

% Kite forces
W_K = p.m_K*p.g*k_E; % in S_E
W_K = R.R_KE*W_K; % in S_K
T_AKp = -T_ASp'; % in S_E
T_AKp = R.R_KE*T_AKp; % in S_K
T_AKm = - T_ASm'; % in S_E
T_AKm = R.R_KE*T_AKm; % in S_K
F_K = W_K + T_AKp + T_AKm + F_a; % in S_K

% Kite moment about OK
M_OK = cross((R.R_KE*OK_AKp'),T_AKp) + cross((R.R_KE*OK_AKm'),T_AKm) + M_a;  % in S_K
end