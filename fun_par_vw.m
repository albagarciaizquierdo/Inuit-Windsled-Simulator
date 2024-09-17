function [x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_vw(p,points,x_red0_design,frame,vw_max)
%% Description: 
% This function gives the state vector, wind speed, angles of elevation and 
% attack, tether's length and minimum wind speed for values of wind 
% speed from 0 to vw_max 
%% Inputs: 
% p --> struct containing parameters
% points --> number of points used in "fun_par_lpoints"
% x_red0_design --> initial guess for the design point
% frame --> reference frame to show kite's results 
%   frame = 1 --> SK
%   frame = 0 --> SE
% vw_max --> maximum wind speed to be calculated
%% Ouputs:
% x_eq --> equilibrium state matrix 
% v_w --> wind speed vector
% gamma --> elevation angle vector
% alpha --> angle of attack vector
% l --> tether's length vector
% v_w_min --> minimum wind speed (scalar)

%%

[x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_lpoints(p,points,x_red0_design,frame);
[x_eq,v_w,gamma,alpha,l] = fun_par_mpoints(p,points,x_red0_design,frame,vw_max,x_eq,v_w,gamma,alpha,l);

end