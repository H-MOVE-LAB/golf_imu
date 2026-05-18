


% Definire la funzione di vincolo
function [c, ceq] = puttingConstraints(par, acc1, acc2, gyr1, gyr2, initial_static, final_static)

c = [];
ceq = [];

intStatic_tot = [initial_static final_static];

beta1 = par(1); 
beta2 = par(2);
bias_acc_1 = par(9:11);
bias_acc_2 = par(12:14);
bias_gyr_1 = par(3:5);
bias_gyr_2 = par(6:8);


acc_1_norm = norm(mean(acc1(intStatic_tot,:) - bias_acc_1));
acc_2_norm = norm(mean(acc2(intStatic_tot,:) - bias_acc_2));

%% Orientation, velocity, and position estimation
intStatic = 1:10;
fs = 100;
N= length(acc1); 
% compute algebraic initial quaternion - initQuaternion(ACC (1 x 3) - bias_acc, MAG (1 x 3))
[qin_1, ~] = initQuaternion( mean(acc1(intStatic,:)) - bias_acc_1,  [0 1 0] );
[qin_2, ~] = initQuaternion( mean(acc2(intStatic,:)) - bias_acc_2,  [1 0 0] );
% class constructor
q1or = MadgwickAHRS('SamplePeriod', 1/fs,'Beta',beta1);
q2or = MadgwickAHRS('SamplePeriod', 1/fs,'Beta',beta2);
% matrix initialization
q_imu1_G = zeros(N,4);
q_imu2_G = zeros(N,4);
q1or.Quaternion = quatconj(qin_1);
q2or.Quaternion = quatconj(qin_2);
q_imu1_G(1,:) = q1or.Quaternion; 
q_imu2_G(1,:) = q2or.Quaternion; 
for jj = 1:N
    % Update 
    q1or.UpdateIMU( gyr1(jj,:) - bias_gyr_1, acc1(jj, :) - bias_acc_1);
    q2or.UpdateIMU( gyr2(jj,:) - bias_gyr_2, acc2(jj, :) - bias_acc_2);
    % Get the updated orientation
    q_imu1_G(jj, :) = q1or.Quaternion;
    q_imu2_G(jj, :) = q2or.Quaternion;
end
% Quaternions -> Rotation matrix
Rimu1_G = quatrotmatr(q_imu1_G);
Rimu2_G = quatrotmatr(q_imu2_G);

t = (0:1/fs:length(acc1)/fs-1/fs)';
g_global = [0; 0; -9.81];

% HEAD
a_measured = acc2 - bias_acc_2; 
a_g_free = zeros(N , 3);
for i = 1:N
    R = Rimu2_G(:,:,i);
    a_g_free(i , :) = (R * a_measured(i , :)' + g_global)';
end
a_2 = a_g_free;
v_2 = cumtrapz(t,a_2);
p_2 = cumtrapz(t,v_2);


% SHAFT
a_measured = acc1 - bias_acc_1; 
a_g_free = zeros(N , 3);
for i = 1:N
    R = Rimu1_G(:,:,i);
    a_g_free(i , :) = (R * a_measured(i , :)' + g_global)';
end
a_1 = a_g_free;
v_1 = cumtrapz(t,a_1);
% p_1 = cumtrapz(t,v_1);



%% Acceleration constraints
% 9.80 < acc[m/s^2] < 9.82 

a = -acc_1_norm + 9.80 ;
c = [c,a];
a =  acc_1_norm - 9.82 ;
c = [c,a];
a = -acc_2_norm + 9.80;
c = [c,a];
a =  acc_2_norm - 9.82;
c = [c,a];


% mean(vel)[m/s] = 0
v(1:3) = mean(v_1(intStatic_tot,:));
ceq = [ceq,v];

v(1:3) = mean(v_2(intStatic_tot,:));
ceq = [ceq,v];

% mean(acc[m/s^2]) = 0

a(1:3) = mean(a_1(intStatic_tot,:));
ceq = [ceq,a];

a(1:3) = mean(a_2(intStatic_tot,:));
ceq = [ceq,a];



%% Velocity constraints
% -0.005 < vel_norm[m/s] < 0.005 
vel_1_norm = norm(mean(v_1(intStatic_tot,:)));
vel_2_norm = norm(mean(v_2(intStatic_tot,:)));

v_inf = -vel_1_norm + 0.005 ;
v_sup =  vel_1_norm - 0.005 ;
c = [ c , v_inf , v_sup ];

v_inf = -vel_2_norm + 0.005 ;
v_sup =  vel_2_norm - 0.005 ;
c = [ c , v_inf , v_sup ];


% -0.005 < vel_max[m/s] < 0.005 
vel_1_max = max(abs(v_1(intStatic_tot,:)));
vel_2_max = max(abs(v_2(intStatic_tot,:)));

v(1:3) = vel_1_max - 0.005 ;
c = [c,v];
v(1:3) = vel_2_max - 0.005 ;
c = [c,v];


% % mean(vel)[m/s] = 0
% v(1:3) = mean(v_1(intStatic_tot,:));
% ceq = [ceq,v];
% 
% v(1:3) = mean(v_2(intStatic_tot,:));
% ceq = [ceq,v];



%% Position constraint
% Head final vertical position < 10 cm 
pz_sup = max(p_2(final_static,3)) - 0.100;
pz_inf = - min(p_2(:,3)) + 0.000;
c = [c,pz_inf,pz_sup];

end