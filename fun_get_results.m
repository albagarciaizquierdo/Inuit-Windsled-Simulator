function [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m,dxdt] = fun_get_results(x,p,frame)
%% Description: 
% This function gives all the results (forces, moments, velocities, 
% positions, lengths and derivative of the state vector) for a given state 
% vector or matrix x
%% Inputs: 
% x --> state vector/matrix
% p --> struct containing parameters
% frame --> reference frame to show kite's results 
%   frame = 1 --> SK
%   frame = 0 --> SE
%% Outputs: 
% F_S --> Sled force
% M_OS --> Sled moment about OS
% T_ASp --> Tether tension from AS+
% T_ASm --> Tether tension from AS-
% W_S --> Sled weight
% N --> Normal force
% F_r --> Friction force
% v_OS --> Sled velocity
% ASp_AKp --> Position vector from AS+ to AK+
% ASm_AKm --> Position vector from AS- to AK-
% up --> unitary vector with origin at AS+ and tips at AK+
% um --> unitary vector with origin at AS- and tips at AK-
% ASp_OS --> Position vector from AS+ to OS
% ASm_OS --> Position vector from AS- to OS
% OK_AKp --> Position vector from OK to AK+
% OK_AKm --> Position vector from OK to AK-

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

% l_p --> tether + length
% l_m --> tether - length

% dxdt --> derivative of the state vector (RHS)
%%
% x should be a matrix nx18
size_x = size(x);
if size_x(2) ~= 18
    x = x';
end

% n --> time steps
n = height(x);

F_S = zeros(n,3);
M_OS = zeros(n,3);
T_ASp = zeros(n,3);
T_ASm = zeros(n,3);
W_S = zeros(n,3);
N = zeros(n,3);
F_r = zeros(n,3);
ASp_AKp = zeros(n,3);
ASm_AKm = zeros(n,3);
up = zeros(n,3);
um = zeros(n,3);
ASp_OS = zeros(n,3);
ASm_OS = zeros(n,3);
OK_AKp = zeros(n,3);
OK_AKm = zeros(n,3);
F_K = zeros(n,3);
M_OK = zeros(n,3);
v_OK = zeros(n,3);
omega_KE = zeros(n,3);
H_OK = zeros(n,3);
W_K = zeros(n,3);
T_AKp = zeros(n,3);
T_AKm = zeros(n,3);
F_a = zeros(n,3);
M_a = zeros(n,3);
alpha = zeros(n,1);
beta = zeros(n,1);
v_a = zeros(n,3);
l_p = zeros(n,1);
l_m = zeros(n,1);
v_OS = zeros(n,3);

dxdt = zeros(18,n);

%% Vectors in S_E
i_E = [1; 0; 0];
j_E = [0; 1; 0];
k_E = [0; 0; 1];
%% Vectors in S_K
i_K = [1 0 0];
j_K = [0 1 0];
k_K = [0 0 1];

for i = 1:n
    % Rotation matrix for the loop
    [R] = fun_rot(x(i,:));
    
    % Sled and kite variables
    [F_S(i,:),M_OS(i,:),T_ASp(i,:),T_ASm(i,:),W_S(i,:),N(i,:),F_r(i,:),v_OS(i,:),ASp_AKp(i,:),ASm_AKm(i,:),up(i,:),um(i,:),ASp_OS(i,:),ASm_OS(i,:),OK_AKp(i,:),OK_AKm(i,:)] = fun_sled(p,x(i,:),R);
    [F_K(i,:),M_OK(i,:),v_OK(i,:),omega_KE(i,:),H_OK(i,:),W_K(i,:),T_AKp(i,:),T_AKm(i,:),F_a(i,:),M_a(i,:),alpha(i,:),beta(i,:),v_a(i,:)] = fun_kite(p,x(i,:),R,OK_AKp(i,:),OK_AKm(i,:),T_ASp(i,:),T_ASm(i,:));
    l_p(i,:) = norm(ASp_AKp(i,:));
    l_m(i,:) = norm(ASm_AKm(i,:));
    
    % Kite variables transformation to SE in case frame ==0
    if frame == 0 
    v_a(i,:) = R.R_EK*v_a(i,:)';
    v_a(i,:) = v_a(i,:)';
    F_K(i,:) = R.R_EK*F_K(i,:)';
    F_K(i,:) = F_K(i,:)';
    F_a(i,:) = R.R_EK*F_a(i,:)';
    F_a(i,:) = F_a(i,:)';
    W_K(i,:) = R.R_EK*W_K(i,:)';
    W_K(i,:) = W_K(i,:)';
    T_AKp(i,:) = R.R_EK*T_AKp(i,:)';
    T_AKp(i,:) = T_AKp(i,:)';
    T_AKm(i,:) = R.R_EK*T_AKm(i,:)';
    T_AKm(i,:) = T_AKm(i,:)';
    M_OK(i,:) = R.R_EK*M_OK(i,:)';
    M_OK(i,:) = M_OK(i,:)';
    v_OK(i,:) = R.R_EK*v_OK(i,:)';
    v_OK(i,:) = v_OK(i,:)';
    omega_KE(i,:) = R.R_EK*omega_KE(i,:)';
    omega_KE(i,:) = omega_KE(i,:)';
    H_OK(i,:) = R.R_EK*H_OK(i,:)';
    H_OK(i,:) = H_OK(i,:)';
    end
        
    % RHS calculation
    dxdt(1,i) = x(i,4);
    dxdt(2,i) = x(i,5);
    dxdt(3,i) = x(i,6);
    dxdt(4,i) = dot(F_S(i,:)'/p.m_S,i_E);
    dxdt(5,i) = dot(F_S(i,:)'/p.m_S,j_E);
    dxdt(6,i) = dot(M_OS(i,:)'/p.Iz_S,k_E);
    dxdt(7,i) = dot((R.R_EK*v_OK(i,:)'),i_E);
    dxdt(8,i) = dot((R.R_EK*v_OK(i,:)'),j_E);
    dxdt(9,i) = dot((R.R_EK*v_OK(i,:)'),k_E);
    dxdt(10,i) = x(i,16) + (x(i,17)*sin(x(i,10))+x(i,18)*cos(x(i,10)))*tan(x(i,11));
    dxdt(11,i) = x(i,17)*cos(x(i,10))-x(i,18)*sin(x(i,10));
    dxdt(12,i) = (x(i,17)*sin(x(i,10))+x(i,18)*cos(x(i,10)))*sec(x(i,11));
    dxdt(13,i) = dot((F_K(i,:)/p.m_K - cross(omega_KE(i,:),v_OK(i,:))),i_K);
    dxdt(14,i) = dot((F_K(i,:)/p.m_K - cross(omega_KE(i,:),v_OK(i,:))),j_K);
    dxdt(15,i) = dot((F_K(i,:)/p.m_K - cross(omega_KE(i,:),v_OK(i,:))),k_K);
    dxdt(16,i) = dot((M_OK(i,:)/p.Ix_K - cross(omega_KE(i,:),H_OK(i,:))),i_K);
    dxdt(17,i) = dot((M_OK(i,:)/p.Iy_K - cross(omega_KE(i,:),H_OK(i,:))),j_K);
    dxdt(18,i) = dot((M_OK(i,:)/p.Iz_K - cross(omega_KE(i,:),H_OK(i,:))),k_K);
    
end

end
