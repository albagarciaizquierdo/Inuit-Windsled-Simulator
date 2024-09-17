function [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm] = fun_sled(p,x,R)
%% Description: 
% This function computes sled forces, moments, velocities and positions in S_E
%% Inputs: 
% x --> state vector/matrix
% p --> struct containing parameters
% R --> struct containing rotation matrices
%% Outputs: (all in S_E in SI units)
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

%% Vectors in S_E
i_E = [1; 0; 0];
j_E = [0; 1; 0];
k_E = [0; 0; 1];

%% Input from tethers
[ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,T_ASp,T_ASm] = fun_tethers(p,x,R);

%% Sled velocity
v_OS = x(4)*i_E + x(5)*j_E;

%% Sled forces
% Weight
W_S = p.m_S*p.g*k_E; % in S_E
% Normal force
N = dot(-W_S - T_ASp - T_ASm,k_E)*k_E; % Normal force
% Friction force
F_r_smax = p.mu_s*norm(N)*i_E;
if norm(v_OS) ==0  % Static friction
    F_r = dot(-T_ASp -T_ASm,i_E)*i_E; % If sled is not moving friction force will oppose potential direction of movement (i.e. wind velocity direction: -i_E)
else % Kinetic friction
    F_r = -p.mu_k*norm(N)*v_OS/norm(v_OS);
    if abs(dot(F_r,v_OS/norm(v_OS))) >= abs(dot(T_ASp,v_OS/norm(v_OS))+dot(T_ASm,v_OS/norm(v_OS)))
        x(4) = 0;
        x(5) = 0;
        v_OS = x(4)*i_E + x(5)*j_E;
        F_r = dot(-T_ASp -T_ASm,i_E)*i_E;
    end

end

% SLED forces
F_S = W_S + N + F_r + T_ASp + T_ASm; % in S_E frame

%% Moment about OS
OS_ASp = -ASp_OS;
OS_ASm = -ASm_OS;
M_OS = cross(OS_ASp,T_ASp) + cross(OS_ASm,T_ASm); % in S_E frame
end