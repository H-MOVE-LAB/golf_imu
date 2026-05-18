function F = minimizeOrientation(par, acc1, acc2, gyr1, gyr2, R_imu2_imu1)

%% Parametri da ottimizzare
beta_1 = par(1);
beta_2 = par(2);
bias_1 = par(3:5);
bias_2 = par(6:8);


%% Madgwick
% orientation computation
intStatic = 1:10;
fs = 100;
N= size(acc1,1); 

% compute algebraic initial quaternion
[qin_1, ~] = initQuaternion( mean(acc1(intStatic,:)) ,  [0 1 0] );
[qin_2, ~] = initQuaternion( mean(acc2(intStatic,:)) ,  [1 0 0] );

% class constructor
q1or = MadgwickAHRS('SamplePeriod', 1/fs,'Beta',beta_1);
q2or = MadgwickAHRS('SamplePeriod', 1/fs,'Beta',beta_2);

% matrix initialization
q_imu1_G = zeros(N,4);
q_imu2_G = zeros(N,4);

q1or.Quaternion = quatconj(qin_1);
q2or.Quaternion = quatconj(qin_2);

q_imu1_G(1,:) = q1or.Quaternion; 
q_imu2_G(1,:) = q2or.Quaternion; 


for jj = 1:N
    % Update 
    q1or.UpdateIMU( gyr1(jj,:) - bias_1, acc1(jj, :) );
    q2or.UpdateIMU( gyr2(jj,:) - bias_2, acc2(jj, :) );
    % Get the updated orientation
    q_imu1_G(jj, :) = q1or.Quaternion;
    q_imu2_G(jj, :) = q2or.Quaternion;
end

% Quaternions -> Rotation matrix
Rimu1_G = quatrotmatr(q_imu1_G);
Rimu2_G = quatrotmatr(q_imu2_G);

% shaft orientation virtually aligned to head and described with respect to the
% global reference
R_imu2_imu1_G= multiprod(Rimu1_G, R_imu2_imu1, [1 2], [1 2]);

q_imu2_G = quatconj(dcm2quat(Rimu2_G));
q_imu2_imu1_G = quatconj(dcm2quat(R_imu2_imu1_G));

q_difference = quatmultiply(quatconj(q_imu2_G), q_imu2_imu1_G);

ang_difference = squeeze(dcm2axangFAST(quatrotmatr(q_difference),false));

ang_difference = rad2deg(ang_difference);



%% Objective function
F = rms(ang_difference)^2;

end
