function [x_eq,v_w,gamma,alpha,l] = fun_par_mpoints(p,points,x_red0_design,frame,vw_max,x_eq,v_w,gamma,alpha,l)
%% Description: 
% This function complements "fun_par_lpoints" giving the state vector, wind 
% speed, angles of elevation and attack, tether's length for values of wind 
% speed larger than the design point p.vw
%% Inputs: 
% p --> struct containing parameters
% points --> number of points used in "fun_par_lpoints"
% x_red0_design --> initial guess for the design point
% frame --> reference frame to show kite's results 
%   frame = 1 --> SK
%   frame = 0 --> SE
% vw_max --> maximum wind speed to be calculated
% x_eq,v_w,gamma,alpha,l --> ouput vector and matrices from "fun_par_lpoints"
%% Ouputs:
% x_eq --> equilibrium state matrix 
% v_w --> wind speed vector
% gamma --> elevation angle vector
% alpha --> angle of attack vector
% l --> tether's length vector

%%
points2 = 1.4*points;
global p
p.v_w = 10;

for i=(length(v_w)):points2
    if i == (length(v_w))
        x_red0(:,i) = x_red0_design;
        v_w(i) = p.v_w;
    else
        x_red0(:,i) = X_red_eq(:,i-1);
        p.v_w = p.v_w + (vw_max-10)/(points2-points);
        v_w(i) = p.v_w;
    end
    [X_red_eq(:,i),~,~] = my_fzero("fun_equilibrio_red",x_red0(:,i),1e-8,30,1e-6);
    x_eq(:,i) = fun_equilibrium_conditions(X_red_eq(:,i));
    [gamma(i), ~] = fun_gamma(x_eq(:,i),p);
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,alpha(i),~,~,l(:,i),~,~] = fun_get_results(x_eq(:,i),p,frame);
end