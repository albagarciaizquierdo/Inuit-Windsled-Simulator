function R_EK = fun_R_EK(phi,theta,psi)
% Description: This function get the rotation matrix that transforms
% vectors from frame S_K to S_E
% Input: angles psi (yaw angle), theta (pitch angle), phi (roll angle)  in radians
% Output: rotation matrix R_EK

R_EK = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
        cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
        -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];

end