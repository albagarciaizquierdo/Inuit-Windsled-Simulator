function g  = fun_equilibrio_paper(t,alpha_eq)
%% Description: 
% This function gives the function g from G. Sánchez Arriaga, M. García
% Villalba and R. Schmehl "Modeling and dynamics of a two-line kite" for
% the equilibrium angle of attack
%% Inputs: 
% t --> time
% alpha_eq --> equilibrium angle of attack
%% Outputs: 
% g --> function 

%%
global p
g = p.Lambda*(p.xA_K*cos(alpha_eq)+p.zA_K*sin(alpha_eq)) + p.c*(p.Cm0+p.Cmalpha*alpha_eq) - p.zA_K*(p.Cx0+p.Cxalpha*alpha_eq) + p.xA_K*(p.Cz0+p.Czalpha*alpha_eq);

end