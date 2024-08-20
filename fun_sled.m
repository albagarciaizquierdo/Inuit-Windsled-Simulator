function [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,v_OS,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm] = fun_sled(p,x,R)
%% Description: This function computes forces and moments acting ond the
% sled
% Inputs: parameters, state vector, output from fun_tethers
% Outputs: Sled forces and moments in S_E

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
    fprintf('Static part \n')
    F_r = dot(-T_ASp -T_ASm,i_E)*i_E; % If sled is not moving friction force will oppose potential direction of movement (i.e. wind velocity direction: -i_E)
else % Kinetic friction
%     fprintf('Kinetic part ')
    F_r = -p.mu_k*norm(N)*v_OS/norm(v_OS);
    if abs(dot(F_r,v_OS/norm(v_OS))) >= abs(dot(T_ASp,v_OS/norm(v_OS))+dot(T_ASm,v_OS/norm(v_OS)))
        fprintf('Kinetic but Fr>T \n')
        x(4) = 0;
        x(5) = 0;
        v_OS = x(4)*i_E + x(5)*j_E;
        F_r = dot(-T_ASp -T_ASm,i_E)*i_E;
    else
        fprintf('Kinetic still \n')
    end

end
% fprintf('Fr sled = ');
% fprintf('%.2f ',F_r);
% fprintf('\n')

% SLED forces
F_S = W_S + N + F_r + T_ASp + T_ASm; % in S_E frame

%% Moment about OS
OS_ASp = -ASp_OS;
OS_ASm = -ASm_OS;
M_OS = cross(OS_ASp,T_ASp) + cross(OS_ASm,T_ASm); % in S_E frame
end