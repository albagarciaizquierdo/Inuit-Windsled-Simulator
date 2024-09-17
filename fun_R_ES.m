 function R_ES = fun_R_ES(psi_S)
%% Description: 
% This function gets the rotation matrix that transforms
% vectors from frame S_S to S_E
%% Input: 
% psi_S --> yaw angle in radians
%% Output: 
% R_ES --> rotation matrix that transforms vectors from frame S_S to S_E

%%
R_ES = [cos(psi_S) -sin(psi_S) 0;
        sin(psi_S) cos(psi_S)  0;
        0   0   1];

end