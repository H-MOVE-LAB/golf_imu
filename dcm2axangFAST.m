function [ang, ax] = dcm2axangFAST(R, flagAxis)
% Function to fast compute the axis-angle representation statring from the
% corresponding rotation matrix.

% INPUT
% R = [3x3] rotation matrix
% flagAxis = logical flag, set to FALSE to output only the angle representation or to TRUE to output both the angle and axis representation

% OUTPUT
% ang = angular value (radians)
% ax = [3x1] unit vector representing the rotation axis

% B. Siciliano, L. Sciavicco, L. Villani, and G. Oriolo, “Robotics,” 2009, doi: 10.1007/978-1-84628-642-1.
% Implemented by Marco Caruso (marco.caruso@polito.it) and Elisa Digo (elisa.digo@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy;
% Department of Mechanical and Aerospace Engineering, Politecnico di Torino, Torino, Italy;
% Last modified: 10/04/2024

r11 = R(1,1,:);
r22 = R(2,2,:);
r33 = R(3,3,:);

ang = acos((r11 + r22 + r33 - 1)/2);

if flagAxis
    r32 = R(3,2,:);
    r23 = R(2,3,:);
    r13 = R(1,3,:);
    r31 = R(3,1,:);
    r21 = R(2,1,:);
    r12 = R(1,2,:);

    ax = 1/(2*sin(ang))*[r32 - r23; r13 - r31; r21 - r12];
end
end
