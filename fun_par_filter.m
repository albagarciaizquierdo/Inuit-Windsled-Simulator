function [v_w_filt,gamma_filt,alpha_filt,l_filt,x_eq_filt] = fun_par_filter(v_w_min,v_w,gamma,alpha,l,x_eq)
%% Description: 
% This function filters the wind speed, angles of elevation and attack, 
% length and equilibrium state vector for the physical values (i.e. when
% v_w > v_w_min)
%% Inputs: 
% v_w_min --> minimum wind speed (scalar)
% v_w --> wind speed vector
% gamma --> elevation angle
% alpha --> angle of attack 
% l --> tether's length
% x_eq --> equilibrium state vector
%% Outputs:
% filtered inputs

%%
% Find indices for which v_w is bigger than v_w_min
valid_indices = v_w >= v_w_min;

% Filter 
v_w_filt = v_w(valid_indices);
gamma_filt = gamma(valid_indices);
alpha_filt = alpha(valid_indices);
l_filt = l(valid_indices);
x_eq_filt = x_eq(:,valid_indices);

end