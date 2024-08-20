function [R] = fun_rot(x)
% Description: This function collects all the necessary rotation matrices
% Input: state vector "x"
% Ouput: Struct "R" with rotation matrices R_ES, R_EK, R_KE

R.R_ES = fun_R_ES(x(3)); % Converts vectors from S_S to S_E
R.R_EK = fun_R_EK(x(10),x(11),x(12)); % Converts vectors from S_K to S_E
R.R_KE = R.R_EK'; % Converts vectors from S_E to S_K
end