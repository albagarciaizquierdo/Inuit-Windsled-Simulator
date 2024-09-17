function [ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,T_ASp,T_ASm] = fun_tethers(p,x,R)
%% Description: 
% This function computes tether vectors and tensions
%% Inputs: 
% x --> state vector/matrix
% p --> struct containing parameters
% R --> struct containing rotation matrices
%% Outputs: (all in S_E in SI units)
% ASp_AKp --> Position vector from AS+ to AK+
% ASm_AKm --> Position vector from AS- to AK-
% up --> unitary vector with origin at AS+ and tips at AK+
% um --> unitary vector with origin at AS- and tips at AK-
% ASp_OS --> Position vector from AS+ to OS
% ASm_OS --> Position vector from AS- to OS
% OK_AKp --> Position vector from OK to AK+
% OK_AKm --> Position vector from OK to AK-
% T_ASp --> Tether tension from AS+
% T_ASm --> Tether tension from AS-

%% Vectors in S_E
i_E = [1; 0; 0];
j_E = [0; 1; 0];
k_E = [0; 0; 1];

% Position vector of the kite Center of mass
OE_OK = x(7)*i_E + x(8)*j_E + x(9)*k_E;

% Position vector of the sled Center of Mass
OE_OS = x(1)*i_E + x(2)*j_E; 
OS_OE = -OE_OS;
%% Tethers
% Tether +
ASp_OS  = R.R_ES*[-p.xA_S; -p.yA_S; 0]; % in S_E frame
OK_AKp  = R.R_EK*[p.xA_K; p.yA_K; p.zA_K]; % in S_E frame
ASp_AKp = ASp_OS + OS_OE + OE_OK + OK_AKp; % in S_E frame
up      = ASp_AKp/norm(ASp_AKp); % in S_E frame
eps_p = norm(ASp_AKp)/p.l0 - 1;

% Tether -
ASm_OS  = R.R_ES*[-p.xA_S; p.yA_S; 0]; % in S_E frame
OK_AKm  = R.R_EK*[p.xA_K; -p.yA_K; p.zA_K]; % in S_E frame
ASm_AKm = ASm_OS + OS_OE + OE_OK + OK_AKm; % in S_E frame
um      = ASm_AKm/norm(ASm_AKm); % in S_E frame
eps_m = norm(ASm_AKm)/p.l0 - 1;

%% Tensions
if norm(ASp_AKp)<p.l0
    T_ASp = zeros(3,1);
else
    T_ASp = p.k*p.l0*eps_p*up; % Tension from AS+ in S_E frame
end

%% Tensions
if norm(ASm_AKm)<p.l0
    T_ASm = zeros(3,1);
else
    T_ASm = p.k*p.l0*eps_m*um; % Tension from AS- in S_E frame
end


end


