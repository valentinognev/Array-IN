%% Single Trajectory simulation 
%% 
% 
% Simulations
% Create noisy measurements and a complicated trajectory.
clear all
% close all
 
showfigures=false;

pathScript = fileparts(matlab.desktop.editor.getActiveFilename);
pathData = fullfile(pathScript,"data");
mkdir(pathData);
%%
% save(fullfile(pathData,"workspace.mat"))
%%
% Run
Fs = 500;
% Fs = 100;
Ts = 1/Fs;  % [s] Sampling period
time_factor = 4.5;
T_end = 10*time_factor;  % [s]
t = 0:Ts:T_end - Ts;
N = length(t);

F_pos = 100;
assert(mod(Fs,F_pos) == 0)
N_pos_update = Fs/F_pos;
inds_pos_update = 1:N_pos_update:N;

% Low dynamics
% N_periods = 1.5 % around 2000 deg/s in norm(w)
% High dynamics
N_periods = 10; % Around 1000 deg/s in norm(w)

time_point_release = 40;

sig_cont_gyro_deg = 1/sqrt(500); % Gyro continous noise 
sig_cont_gyro = deg2rad(sig_cont_gyro_deg); 
sig_cont_acc = 0.5/sqrt(500);     % Acc continous noise 

sig_acc = sqrt(Fs)*sig_cont_acc;     % Discrete Acc noise
sig_gyro = sqrt(Fs)*sig_cont_gyro;   % Discrete gyro noise
rad2deg(sig_gyro)

sig_pos = 1e-1;        % Sigma position update

% estimate a constant bias
sig_b_a = 0;
sig_b_g = 0;
sig_b_a_init = 0.2;          % accelerometer bias value, From  Carlsson2021
sig_b_g_init = deg2rad(1);   % gyro bias value, From  Carlsson2021


% phi in [0, pi]
inp.phi = "sinus";
params = struct;
params.A = pi/2;
params.b = pi/2;
params.f = (2 + N_periods)/T_end*time_factor; % Three periods
inp.phi_params = params;

% theta in [0, 2pi]
inp.theta = "sinus";
params = struct;
params.A = pi;
params.b = pi;
params.f = (N_periods)/T_end*time_factor; % two periods
inp.theta_params = params;

m = get_spherical_motion(t, inp);   % calculate phi,theta, phi_dot, theta_dot history

vels = compute_omega(m); % calculation of omega and omega_dot in navigation and body frames

% f = figure(1);
% set(f,'WindowStyle','docked')
% clf

if showfigures
    figure(1);
    vels = mergeStruct(m, vels);
    plot_angular_velocities(t, vels)
    scale_figure(gcf, 1.5);
    figure
    hold on
    plot(t, rad2deg(norm_time(vels.w_b)), "b-")
    plot(t, rad2deg(norm_time(vels.w_n)), "r--")
end
%%

phi = m.phi;
theta = m.theta;


R_bn_true = R_bn_spherical(phi,theta); % calculate navigation to body frame rotation matrix history
R_offset = eye(3);

w_b = R_offset*vels.w_b;
w_dot_b = R_offset*vels.w_dot_b;


% Make initial rotation eye(3)
for n = 1:N
    R_bn_true(:,:,n) = R_offset'*R_bn_true(:,:,n);
end
R_nb_true = zeros(size(R_bn_true));
for n = 1:N
    R_nb_true(:,:,n) = R_bn_true(:,:,n)';
end
%%

if showfigures    
    figure(2);
    eul = my_rotm2eul(R_nb_true);
    plot(t,rad2deg(eul'))
    legend("roll", "pitch", "yaw")
    grid on
end

%% sensor positions in body frame
pos_b = [-0.0095    0.0032    0.0010
   -0.0095    0.0032   -0.0010
   -0.0095    0.0095    0.0010
   -0.0095    0.0095   -0.0010
   -0.0095   -0.0095    0.0010
   -0.0095   -0.0095   -0.0010
   -0.0095   -0.0032    0.0010
   -0.0095   -0.0032   -0.0010
   -0.0032    0.0032    0.0010
   -0.0032    0.0032   -0.0010
   -0.0032    0.0095    0.0010
   -0.0032    0.0095   -0.0010
   -0.0032   -0.0095    0.0010
   -0.0032   -0.0095   -0.0010
   -0.0032   -0.0032    0.0010
   -0.0032   -0.0032   -0.0010
    0.0032    0.0032    0.0010
    0.0032    0.0032   -0.0010
    0.0032    0.0095    0.0010
    0.0032    0.0095   -0.0010
    0.0032   -0.0095    0.0010
    0.0032   -0.0095   -0.0010
    0.0032   -0.0032    0.0010
    0.0032   -0.0032   -0.0010
    0.0095    0.0032    0.0010
    0.0095    0.0032   -0.0010
    0.0095    0.0095    0.0010
    0.0095    0.0095   -0.0010
    0.0095   -0.0095    0.0010
    0.0095   -0.0095   -0.0010
    0.0095   -0.0032    0.0010
    0.0095   -0.0032   -0.0010]';
% i_IMU_select = [3,5,27,29];
% i_IMU_select = sort([i_IMU_select, i_IMU_select+1])
i_IMU_select = 1:32;
pos_b = pos_b(:,i_IMU_select);   % Select specific IMUs for calculations
%%
r_tot = pos_b - mean(pos_b,2);   % Centroid of sensor positions

%%
if showfigures
    figure(3);
    subplot(1,2,1)
    plot(r_tot(1,1:2:end)*1e3,r_tot(2,1:2:end)*1e3,"x")
    axis("equal")
    grid on
    L = 15;
    xlim([-L,L])
    ylim([-L,L])
    % xlim([-5,5]*1e-3)
    xlabel("[mm]")
    ylabel("[mm]")
    title("Topside")
    subplot(1,2,2)
    plot(r_tot(1,2:2:end)*1e3,r_tot(2,2:2:end)*1e3,"x")
    grid on
    axis("equal")
    xlim([-L,L])
    ylim([-L,L])
    xlabel("[mm]")
    ylabel("[mm]")
    title("underside")
end
%%

K = size(pos_b,2); % Number of acc triads
inds = reshape(1:3*K,3,[]);

settings_default = struct;
settings_default.save_full_covariances = false;
settings_default.save_jacobians = false;
settings_default.verbose = true;
settings_default.T = Ts;                            % Sampling period
settings_default.g = [0; 0; -9.81];


settings_default.r = r_tot;


%% Generate the motion 

f_pos = [1 2 3 ]/T_end*2*pi*time_factor; % frequency for each axis


p_true = zeros(3,N);
v_true = zeros(3,N);
a_true = zeros(3,N);

for i = 1:3
    f_i = f_pos(i);
    p_true(i,:) = sin(f_i*t);
    v_true(i,:) = f_i*cos(f_i*t);
    a_true(i,:) = -f_i^2*sin(f_i*t);
end
%%

if showfigures
    figure(4);
    hold on
    for i = 1:3
        plot(t, p_true(i,:))

    end
    grid on
    legend("x","y","z")
    ylabel("[m]")
    xlabel("Time [s]")
    title("True Navigation Position");

    figure(5)
    hold on
    for i = 1:3
        plot(t, v_true(i,:))

    end
    grid on
    legend("x","y","z")
    ylabel("[m/s]")
    xlabel("Time [s]")
    title("True Navigation Velocity");


    figure(6)
    hold on
    for i = 1:3
        plot(t, a_true(i,:))

    end
    grid on
    legend("x","y","z")
    ylabel("[m/s]")
    xlabel("Time [s]")
    title("True Navigation Acceleration");
end
%% 
% Zero acceleration, stay at the same position. Measure gravity in body frame

u_true_acc = zeros(3*K,N);
rot_comp = zeros(3*K,N);
trans_comp = zeros(3,N);
for n = 1:N
    W = skew_sym(w_b(:,n))^2 + skew_sym(w_dot_b(:,n));
    trans_comp(:,n) = R_bn_true(:,:,n)*(a_true(:,n) - settings_default.g);   % rotate body center accelerations 
                                                                             % from inertial frame to body frame
    for k = 1:K
        kk = inds(:,k);
        r_k = r_tot(:,k);         % sensors position 
        rot_comp(kk,n) = W*r_k;   % rotational acceleration component
        u_true_acc(kk,n) = trans_comp(:,n) + rot_comp(kk,n);  % total acceleration measured in body frame
    end
end

% Gyroscopes measure angular velocity in body coordinates
u_true_gyro = w_b;
%%
rng("default");  % set random number generator type

b_g_true = sig_b_g_init*randn(3,1)/sqrt(K);
b_a_true = sig_b_a_init*randn(3*K,1);            

A = compute_A_non_center(r_tot);
b_omega_dot_true = -A(1:3,:)*b_a_true;    % adding measurement bias to all omega_dot values
b_s_true = -A(4:6,:)*b_a_true;

S_true = struct;
S_true.R = R_nb_true;
S_true.w = w_b;
S_true.p = p_true;
S_true.v = v_true;
S_true.omega_dot = w_dot_b;

S_true.b_omega_dot = repmat(b_omega_dot_true,1,N);  % omega_dot bias 
S_true.b_s = repmat(b_s_true,1,N);
S_true.b_g = repmat(b_g_true,1,N);
S_true.s = trans_comp;                              % body linear acceleration in body frame

Q_acc = eye(3*K)*sig_acc^2;
Q_gyro = eye(3)*sig_gyro^2/K;
Q_pos = eye(3)*sig_pos^2;

noise_acc = chol(Q_acc)*randn(3*K,N);
noise_gyro = chol(Q_gyro)*randn(3,N);
noise_pos = nan(3,N);
% Position update every 10 sample 
noise_pos(:,inds_pos_update) = chol(Q_pos)*randn(3, length(inds_pos_update));

% Release at 15 secs
noise_pos(:,t > time_point_release) = nan;


%%

% Initial values 

rng("default");

init = struct;
init.cov.R = eye(3)*deg2rad(1)^2;
init.cov.p = Q_pos*2; 
init.cov.v = eye(3)*0.1^2;
init.cov.omega = Q_gyro*2;

init.mean.R = R_nb_true(:,:,1)*expSO3(chol(init.cov.R)*randn(3,1));
init.mean.p = p_true(:,1) + noise_pos(:,1); 
init.mean.v = v_true(:,1) + chol(init.cov.v)*randn(3,1);
init.mean.omega = u_true_gyro(:,1) + noise_gyro(:,1);

init.mean.b_a = zeros(3*K,1);
init.cov.b_a = sig_b_a_init^2*eye(3*K)*3;

init.mean.b_g = zeros(3,1);
init.cov.b_g = eye(3)*sig_b_g_init^2/K*3;


%% Run Filters

% error("asd")
% Models

simdata = struct;

simdata.accelerometer_array_2nd_order = settings_default;
simdata.accelerometer_array_2nd_order.propagate_position = true;
simdata.accelerometer_array_2nd_order.propagate_velocity = true;
simdata.accelerometer_array_2nd_order.propagate_bias_alpha = true;
simdata.accelerometer_array_2nd_order.propagate_bias_gyro = true;
simdata.accelerometer_array_2nd_order.set_T2_R_zero = false;
simdata.accelerometer_array_2nd_order.get_model = @D_LG_EKF_Array_v4_alpha;
simdata.accelerometer_array_2nd_order.r = r_tot;
simdata.accelerometer_array_2nd_order.Q_acc = sig_acc^2*eye(3*K);
simdata.accelerometer_array_2nd_order.Q_bias_acc = sig_b_a^2*eye(3*K);
simdata.accelerometer_array_2nd_order.R_pos = sig_pos^2*eye(3);
simdata.accelerometer_array_2nd_order.R_gyro = sig_gyro^2*eye(3)/K;
simdata.accelerometer_array_2nd_order.Q_bias_gyro = sig_b_g^2*eye(3)/K;
simdata.accelerometer_array_2nd_order.do_gyro_updates = true;
simdata.accelerometer_array_2nd_order.do_position_updates = true;
simdata.accelerometer_array_2nd_order.do_rotation_updates = false;
simdata.accelerometer_array_2nd_order.label = "2nd order accelerometer array";


simdata.accelerometer_array_1st_order = settings_default;
simdata.accelerometer_array_1st_order.propagate_position = true;
simdata.accelerometer_array_1st_order.propagate_velocity = true;
simdata.accelerometer_array_1st_order.propagate_bias_alpha = true;
simdata.accelerometer_array_1st_order.propagate_bias_gyro = true;
simdata.accelerometer_array_1st_order.set_T2_R_zero = true;
simdata.accelerometer_array_1st_order.get_model = @D_LG_EKF_Array_v4_alpha;
simdata.accelerometer_array_1st_order.r = r_tot;
simdata.accelerometer_array_1st_order.Q_acc = sig_acc^2*eye(3*K);
simdata.accelerometer_array_1st_order.Q_bias_acc = sig_b_a^2*eye(3*K);
simdata.accelerometer_array_1st_order.R_pos = sig_pos^2*eye(3);
simdata.accelerometer_array_1st_order.R_gyro = sig_gyro^2*eye(3);
simdata.accelerometer_array_1st_order.Q_bias_gyro = sig_b_g^2*eye(3)/K;
simdata.accelerometer_array_1st_order.do_gyro_updates = true;
simdata.accelerometer_array_1st_order.do_position_updates = true;
simdata.accelerometer_array_1st_order.do_rotation_updates = false;
simdata.accelerometer_array_1st_order.label = "1st order accelerometer array";

simdata.gyroscope_2nd_order = settings_default;
simdata.gyroscope_2nd_order.propagate_bias_gyro = true;
simdata.gyroscope_2nd_order.input_accelerometers = true;
simdata.gyroscope_2nd_order.propagate_position = true;
simdata.gyroscope_2nd_order.propagate_velocity = true;
simdata.gyroscope_2nd_order.propagate_bias_s = true;
simdata.gyroscope_2nd_order.set_T2_R_zero = false;
simdata.gyroscope_2nd_order.get_model = @D_LG_EKF_Gyro_2nd_v4;
simdata.gyroscope_2nd_order.r = r_tot;
simdata.gyroscope_2nd_order.Q_acc = sig_acc^2*eye(3*K);
simdata.gyroscope_2nd_order.Q_bias_acc = sig_b_a^2*eye(3*K);
simdata.gyroscope_2nd_order.R_pos = sig_pos^2*eye(3);
simdata.gyroscope_2nd_order.Q_gyro = sig_gyro^2*eye(3)/K;
simdata.gyroscope_2nd_order.Q_bias_gyro = sig_b_g^2*eye(3)/K;
simdata.gyroscope_2nd_order.do_gyro_updates = false;
simdata.gyroscope_2nd_order.do_position_updates = true;
simdata.gyroscope_2nd_order.do_rotation_updates = false;
simdata.gyroscope_2nd_order.label = "2nd order gyroscope";

simdata.gyroscope_1st_order = settings_default;
simdata.gyroscope_1st_order.input_accelerometers = true;
simdata.gyroscope_1st_order.propagate_bias_s = true;
simdata.gyroscope_1st_order.propagate_bias_gyro = true;
simdata.gyroscope_1st_order.propagate_position = true;
simdata.gyroscope_1st_order.propagate_velocity = true;
simdata.gyroscope_1st_order.get_model = @D_LG_EKF_Gyro_1st_v4;
simdata.gyroscope_1st_order.N_a = K;
simdata.gyroscope_1st_order.Q_acc = sig_acc^2*eye(3*K);
simdata.gyroscope_1st_order.Q_bias_acc = sig_b_a^2*eye(3*K);
simdata.gyroscope_1st_order.R_pos = sig_pos^2*eye(3);
simdata.gyroscope_1st_order.Q_gyro = sig_gyro^2*eye(3)/K;
simdata.gyroscope_1st_order.Q_bias_gyro = sig_b_g^2*eye(3)/K;
simdata.gyroscope_1st_order.do_gyro_updates = false;
simdata.gyroscope_1st_order.do_position_updates = true;
simdata.gyroscope_1st_order.do_rotation_updates = false;
simdata.gyroscope_1st_order.label = "1st order gyroscope";

% 


sensorData = struct;

sensorData.acc_measurements = u_true_acc + b_a_true + noise_acc;
sensorData.gyro_measurements = u_true_gyro + b_g_true + noise_gyro;
sensorData.position_measurements = p_true + noise_pos;

%% Run filters
% 

resSingle = struct;

run_settings = struct;
run_settings.compute_error = true;
run_settings.verbose = true;
run_settings.save_input = true;
run_settings.save_residuals = true;

resSingle.label = "";
[~,resSingle.accelerometer_array_1st_order] = run_filter(sensorData, init, simdata.accelerometer_array_1st_order, run_settings, S_true);
[~, resSingle.accelerometer_array_2nd_order] = run_filter(sensorData, init, simdata.accelerometer_array_2nd_order, run_settings, S_true);
[~, resSingle.gyroscope_2nd_order] = run_filter(sensorData, init, simdata.gyroscope_2nd_order,run_settings, S_true);
[~,resSingle.gyroscope_1st_order] = run_filter(sensorData, init, simdata.gyroscope_1st_order,run_settings, S_true);

%%
true_traj = struct;
true_traj.filt.mean.b_g = S_true.b_g;
true_traj.filt.mean.b_omega_dot = S_true.b_omega_dot;
true_traj.filt.mean.b_s = S_true.b_s;
true_traj.filt.mean.p = p_true;
true_traj.filt.mean.v = v_true;
true_traj.filt.mean.R = R_nb_true;
true_traj.filt.mean.w = w_b;
true_traj.filt.std = struct;
true_traj.label = "True";

resSingle.True = struct;
resSingle.True.res = true_traj;
%%
cases_plot = ["mean_array", "mean_array_omega_dot", "array_1st", "array_2nd","True"];

plotOpt = struct;
% plotOpt.mean_array.cov.ls = "-";
plotOpt.mean_array.cov.m = "+";
plotOpt.mean_array.cov.m_space = 1000;
plotOpt.mean_array.cov.m_offset = 1;

plotOpt.mean_array_omega_dot.cov.m = "+";
plotOpt.mean_array_omega_dot.cov.m_space = 1000;
plotOpt.mean_array_omega_dot.cov.m_offset = 100;

plotOpt.array_1st.cov.m = "+";
plotOpt.array_1st.cov.m_space = 1000;
plotOpt.array_1st.cov.m_offset = 200;

plotOpt.array_2nd.cov.m = "+";
plotOpt.array_2nd.cov.m_space = 1000;
plotOpt.array_2nd.cov.m_offset = 300;

% Bias specific force

cases_plot = [ "accelerometer_array_2nd_order";
     "accelerometer_array_1st_order";
     "gyroscope_2nd_order";
     "gyroscope_1st_order";
"True"];


plotOpt_k = struct;
extras = struct;
extras.path = ["res", "filt"];
[fig,~] = plot_bias_s(resSingle, t, cases_plot, plotOpt_k, "normal", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
     "accelerometer_array_1st_order";
     "gyroscope_2nd_order";
     "gyroscope_1st_order"];


plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
[fig,~] = plot_bias_s(resSingle, t, cases_plot, plotOpt_k, "error", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
     "accelerometer_array_1st_order";
     "gyroscope_2nd_order";
     "gyroscope_1st_order"];


plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
extras.show_mean = false;
[fig,~] = plot_bias_s(resSingle, t, cases_plot, plotOpt_k, "error_log", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order"; "True"];


plotOpt_k = struct;
extras = struct;
extras.path = ["res", "filt"];
[fig,~] = plot_bias_omega_dot(resSingle, t, cases_plot, plotOpt_k, "normal", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
     "accelerometer_array_1st_order";
     "gyroscope_2nd_order";
     "gyroscope_1st_order"];


plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
extras.show_mean = false;
[fig,~] = plot_bias_omega_dot(resSingle, t, cases_plot, plotOpt_k, "error_log", extras);
scale_figure(fig,1.5);
% Bias gyroscope

cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order"; "True"];
plotOpt_k = struct;
extras = struct;
extras.path = ["res", "filt"];
[fig,~] = plot_bias_gyroscopes(resSingle, t, cases_plot, plotOpt_k, "normal", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order"];
plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
extras.show_mean = false;
[fig,~] = plot_bias_gyroscopes(resSingle, t, cases_plot, plotOpt_k, "error_log", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order"];
plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
[fig,~] = plot_rotation(resSingle, t, cases_plot, plotOpt_k, "error", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order"; "True"];

plotOpt_k = struct;
extras = struct;
extras.path = ["res", "filt"];
[fig,~] = plot_navigation_position(resSingle, t, cases_plot, plotOpt_k, "normal", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order";];

plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
[fig,~] = plot_navigation_position(resSingle, t, cases_plot, plotOpt_k, "error", extras);
scale_figure(fig,1.5);
%%
cases_plot = [ "accelerometer_array_2nd_order";
    "accelerometer_array_1st_order";
   "gyroscope_2nd_order";
  "gyroscope_1st_order";];

plotOpt_k = struct;
extras = struct;
extras.path = ["res", "err"];
extras.show_mean = false;
[fig,~] = plot_navigation_position(resSingle, t, cases_plot, plotOpt_k, "error_log", extras);
scale_figure(fig,1.5);
%%

%% Functions

function [c, misc] = run_filter(sensorData, init, simdata, run_settings, S_ref, x, f_change, mask_logLL)

if nargin < 4
    run_settings = struct;
end

if isfield(run_settings, "verbose")
    simdata.verbose = run_settings.verbose;
end
if nargin > 5
    [sensorData, init, simdata, changes] = f_change(x, sensorData, init, simdata);
end
model = simdata.get_model(simdata);

res = DLGEKFv4(sensorData, init, model, simdata);

if isfield(run_settings, "compute_error") && nargin > 4
    res.err = compute_error_tot(res, S_ref);
end

misc = struct;
misc.res = res;
if isfield(run_settings, "save_input") && run_settings.save_input
    misc.sensorData = sensorData;
    misc.init = init;
    misc.simdata = simdata;
    misc.model = model;
end

if nargin > 5
    misc.changes = changes;
end

o = 0; % Offset in ind
% Depends on in the order of measurement updates in filter
save_residuals = isfield(run_settings, "save_residuals") && run_settings.save_residuals;

if isfield(sensorData, "gyro_measurements") && simdata.do_gyro_updates && save_residuals 
    mask = all(~isnan(sensorData.gyro_measurements));
    e = res.logL.residuals_normalized(mask);
    e_gyro = cellfun(@(e_i) e_i(1:3), e(2:end), 'UniformOutput', false); % Skip initial value
    misc.e_gyro = cat(2, e_gyro{:});
    o = o + 3;
end

if isfield(sensorData, "position_measurements") && simdata.do_position_updates && save_residuals
    inds = (1:3) + o;
    mask = all(~isnan(sensorData.position_measurements));
    e = misc.res.logL.residuals_normalized(mask);
    e_pos = cellfun(@(e_i) e_i(inds), e(2:end), 'UniformOutput', false); % Skip initial value
    misc.e_pos = cat(2, e_pos{:});
    o = o + 3;
end

if isfield(sensorData, "rotation_measurements") && simdata.do_rotation_updates && save_residuals
    inds = (1:3) + o;
    mask = reshape(all(~isnan(sensorData.rotation_measurements),[1 2]),[],1);
    e = misc.res.logL.residuals_normalized(mask);
    e_rot = cellfun(@(e_i) e_i(inds), e(2:end), 'UniformOutput', false); % Skip initial value
    misc.e_rot = cat(2, e_rot{:});
end

use_mask_logLL = isfield(run_settings, "use_mask_logLL") && run_settings.use_mask_logLL;
if nargin > 7 && use_mask_logLL
    % Assume first time point is included
    parts = res.logL.parts(mask_logLL);
    % Skip first time-point
    logL = -1/2*sum(parts(2:end));
    c = -logL;
else
    c = -res.logL.value;
end

end
%% 
% 

function err = compute_error_tot(S_tot, S_ref)

err = struct;
err.mean = compute_error(S_tot.filt.mean, S_ref);

% Copy standard deviations
fields = fieldnames(err.mean);
for idx = 1:length(fields)
    err.std.(fields{idx}) = S_tot.filt.std.(fields{idx});
end

end
%%