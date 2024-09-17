function dxdt = RHS(t,x)
%% Description: 
% This function gives the Right Hand Side of the derivative of the state 
% vector
%% Input: 
% t --> time (s)
% x --> State vector:
%   x(1)  = x_s,  x(2)  = y_s,     x(3)  = psi_s
%   x(4)  = vx_s, x(5)  = vy_s,    x(6)  = r_s 
%   x(7)  = x_k,  x(8)  = y_k,     x(9)  = z_k
%   x(10) = phi,  x(11) = theta,   x(12) = psi
%   x(11) = u_k,  x(12) = v_k,     x(13) = w_k
%   x(14) = p_k,  x(15) = q_k,     x(16) = r_k
%% Output
% dxdt --> derivative of the state vector

%%
global p  
global frame
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,dxdt] = fun_get_results(x,p,frame);

end

