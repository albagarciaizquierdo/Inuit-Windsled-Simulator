function [v_w_filt,gamma_filt,alpha_filt,l_filt,x_eq_filt] = fun_par_filter(v_w_min,v_w,gamma,alpha,l,x_eq)

% Find indices for which v_w is bigger than v_w_min
valid_indices = v_w >= v_w_min;

% Filter 
v_w_filt = v_w(valid_indices);
gamma_filt = gamma(valid_indices);
alpha_filt = alpha(valid_indices);
l_filt = l(valid_indices);
x_eq_filt = x_eq(:,valid_indices);

end