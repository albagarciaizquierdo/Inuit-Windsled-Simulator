function [gamma_p, gamma_m] = fun_gamma(x,p)
% Gives angle gamma in DEGREES
[R] = fun_rot(x);

[ASp_AKp,ASm_AKm,~,~,~,~,~,~,~,~] = fun_tethers(p,x,R);

k_E = [0;0;1];

gamma_p = asind(dot(ASp_AKp,-k_E)/(norm(ASp_AKp)));
gamma_m = asind(dot(ASm_AKm,-k_E)/(norm(ASm_AKm)));

end