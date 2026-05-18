clearvars
close all
% clc

%% ---------------------------- Data loading ----------------------------
load data_static.mat
load data_trial_1.mat % putting trial

% load putter model
load R_from_shaft_to_head.mat

acc_head = data_dynamic.imu_head.acc;
acc_shaft = data_dynamic.imu_shaft.acc; 
gyr_head = data_dynamic.imu_head.gyr; 
gyr_shaft = data_dynamic.imu_shaft.gyr; 

% Gyro bias (static offset)
offset_gyro_head = mean(data_static.imu_head.gyr); 
offset_gyro_shaft = mean(data_static.imu_shaft.gyr);

%% -----------------------Optimization framework-------------------------
% par_ot = [beta_shaft (1), beta_head (2), delta_Bias_gyr_shaft (3:5),
  % delta_Bias_gyr_head (6:8), delta_Bias_acc_shaft (9:11), delta_Bias_acc_head (12:14)]

% parameter inizialization
par_0 = [0 0 offset_gyro_shaft offset_gyro_head 0 0 0 0 0 0];

% definition of ranges of each parameter
par_lim=sort([
    0   0.2 % beta shaft
    0   0.2 % beta head
    [0.1 1.9].*[offset_gyro_shaft offset_gyro_head]'
    -0.2 0.2 % bias acc shaft
    -0.2 0.2
    -0.2 0.2
    -0.2 0.2 % bias acc head
    -0.2 0.2
    -0.2 0.2
   ],2);

lb = par_lim(:,1);
ub = par_lim(:,2);

initial_static = data_dynamic.initial_static_phase; 
final_static = data_dynamic.final_static_phase; 

optionsFmincon = optimoptions('fmincon','Display','off', 'Algorithm','sqp'); % define the solver options
f_obj =  @(par) minimizeOrientation(par, acc_shaft, acc_head, gyr_shaft, gyr_head, R_head_shaft); % update the objective function
nonlincon = @(par) puttingConstraints(par, acc_shaft, acc_head, gyr_shaft, gyr_head, initial_static, final_static); % non-linear putting-specific constraints
problem = createOptimProblem('fmincon','x0',par_0,'objective',f_obj,'lb',lb,'ub',ub,'options',optionsFmincon,'nonlcon',nonlincon);
ms = MultiStart('Display','off');
[par_opt, res, exflag] = run(ms,problem,7); % find the optimal configuration

%% -------------Standard pipeline using optimized parameters --------------
% Orientation estimation using optimized parameters
intStatic = 1:10;
fs = 100;
N= length(acc_head); 
% compute algebraic initial quaternion - initQuaternion(ACC (1 x 3) - bias_acc , MAG (1 x 3))
[qin_head, ~] = initQuaternion( mean(acc_head(intStatic,:)) - par_opt(12:14) ,  [1 0 0] );
% class constructor
qor_head = MadgwickAHRS('SamplePeriod', 1/fs,'Beta',par_opt(2));
q_head_G = zeros(N,4);
qor_head.Quaternion = quatconj(qin_head);
q_head_G(1,:) = qor_head.Quaternion; 
for jj = 2:N
    % Update - (gyr (1x3) - bias(1x3), acc(1x3) - bias(1x3),  gyroscope units must be radians/second
    qor_head.UpdateIMU( gyr_head(jj,:) - par_opt(6:8), acc_head(jj,:) - par_opt(12:14));
    q_head_G(jj, :) = qor_head.Quaternion;
end
% Quaternions -> Rotation matrix
R_head_G = quatrotmatr(q_head_G);

% Position estimation
g_global = [0; 0; -9.81]; 
a_body_G = zeros(N , 3); 
for jj = 2:N 
    R = R_head_G(:,:,jj); 
    a_body_G(jj , :) = (R * (acc_head(jj,:) - par_opt(12:14))' + g_global)'; 
end

time = 0:1/fs:(N-1)/fs;
v_head = cumtrapz(time,a_body_G);
p_head = cumtrapz(time,v_head);

figure, plot(time,v_head)
figure, plot(time, p_head)

