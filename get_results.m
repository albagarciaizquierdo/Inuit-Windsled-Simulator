function [F_S,M_OS,T_ASp,T_ASm,W_S,N,F_r,ASp_AKp,ASm_AKm,up,um,ASp_OS,ASm_OS,OK_AKp,OK_AKm,F_K,M_OK,v_OK,omega_KE,H_OK,W_K,T_AKp,T_AKm,F_a,M_a,alpha,beta,v_a,l_p,l_m] = get_results(x,p)
% x is matrix nxm
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
T_ASp = zeros(n,3);
T_ASm = zeros(n,3);
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

for i = 1:n
        [F_S(i,:),M_OS(i,:),T_ASp(i,:),T_ASm(i,:),W_S(i,:),N(i,:),F_r(i,:)] = fun_sled(p,x(i,:));
        [ASp_AKp(i,:),ASm_AKm(i,:),up(i,:),um(i,:),ASp_OS(i,:),ASm_OS(i,:),OK_AKp(i,:),OK_AKm(i,:),T_ASp(i,:),T_ASm(i,:)] = fun_tethers(p,x(i,:));
        [F_K(i,:),M_OK(i,:),v_OK(i,:),omega_KE(i,:),H_OK(i,:),W_K(i,:),T_AKp(i,:),T_AKm(i,:),F_a(i,:),M_a(i,:),alpha(i,:),beta(i,:),v_a(i,:)] = fun_kite(p,x(i,:));
        l_p(i,:) = norm(ASp_AKp(i,:));
        l_m(i,:) = norm(ASm_AKm(i,:));
end

end
