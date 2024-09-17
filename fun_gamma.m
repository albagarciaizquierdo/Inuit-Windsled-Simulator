function [gamma_p, gamma_m] = fun_gamma(x,p)
%% Description: 
% This function gives the elevation angles gamma_p, gamma_m in degrees for a given state vector
%% Inputs: 
% x --> state vector (1x18)
% p --> struct containing parameters
%% Outputs: 
% gamma_p, gamma_m --> tethers' elevation angles in degrees

%%
[R] = fun_rot(x);

[ASp_AKp,ASm_AKm,~,~,~,~,~,~,~,~] = fun_tethers(p,x,R);

k_E = [0;0;1];

gamma_p = asind(dot(ASp_AKp,-k_E)/(norm(ASp_AKp)));
gamma_m = asind(dot(ASm_AKm,-k_E)/(norm(ASm_AKm)));

end