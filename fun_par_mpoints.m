function [x_eq,v_w,gamma,alpha,l] = fun_par_mpoints(p,points,x_red0_design,frame,vw_max,x_eq,v_w,gamma,alpha,l)
% function to get More points in the parametric study
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
    x_eq(:,i) = equilibrium_conditions(X_red_eq(:,i));
    [gamma(i), ~] = fun_gamma(x_eq(:,i),p);
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,alpha(i),~,~,l(:,i),~,~] = fun_get_results(x_eq(:,i),p,frame);
end