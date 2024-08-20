function [x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_vw(p,points,x_red0_design,frame,vw_max)

[x_eq, v_w, gamma, alpha, l, v_w_min] = fun_par_lpoints(p,points,x_red0_design,frame);
[x_eq,v_w,gamma,alpha,l] = fun_par_mpoints(p,points,x_red0_design,frame,vw_max,x_eq,v_w,gamma,alpha,l);

end