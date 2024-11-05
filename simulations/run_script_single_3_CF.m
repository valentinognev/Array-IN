%% Single Trajectory simulation
%%
function run_script_single_3_CF
%
% Simulations
% Create noisy measurements and a complicated trajectory.
close all

global showfigures;

showfigures=false;

pathScript = fileparts(matlab.desktop.editor.getActiveFilename);
pathData = fullfile(pathScript,"data");
mkdir(pathData);

[meas, acc, gyro, pos]=initSamplingData();
[phi, theta, vels]=genBodyFrameAngleTrajectory(meas);     % generation of angulat motion
[truth, w_b, w_dot_b] = getNav2BodyRotHistory(phi, theta, meas, vels);
[pos]=getSensorPositions(pos);
[settings_default, inds]=getGeneralSettings(meas, pos);
[truth] = generateMotionTrajectory(meas, truth);          % generation spatial motion
[S_true, acc, gyro, pos, truth] = addNoiseAndBias(meas, pos, settings_default, gyro, acc, w_b, w_dot_b, truth, inds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% START INJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    dataParams.timfactor = 1e-3;    dataParams.posfactor = 1e-3;    dataParams.downsamping = 1;
    data=readtable('../../Data/cf_2004_re00.csv');

    [tim, world, body]=getCFBodyData(data, dataParams);

    
    truth.pos = world.pos;
    truth.v = world.vel;

    ntim=length(tim);
    truth.acc_u=body.acc;
    truth.gyro_u=body.omega;

    meas.N=ntim;
    meas.t=tim;
    meas.T_end = tim(end);
    meas.Ts=tim(2)-tim(1);
    meas.Fs=1/meas.Ts;

    biasFactor = 0.05;   %   1e-8; %  
    acc.init_b_a = biasFactor*max(sqrt(sum(truth.a.^2,1)));
    truth.biasAcc = acc.init_b_a*randn(3,1);

    noiseFactor = 0.01; %  1e-8;   %   
    acc.noiseSigDisc = noiseFactor*max(sqrt(sum(truth.a.^2,1)));
    acc.Q = eye(3)*acc.noiseSigDisc^2;
    acc.noise = chol(acc.Q)*randn(3,meas.N);
    gyro.noise = chol(gyro.Q)*randn(3,meas.N)*0;  % gyro noise is zero by definition, as gyro is absent

    F_pos = meas.Fs/1;                       % [Hz] Time frequency of position update in simulation, must be divider of sampling frequency
    % assert(mod(meas.Fs,F_pos) == 0)
    pos.N_update = floor(meas.Fs/F_pos);     % update position every N steps in simulation
    pos.inds_update = 1:pos.N_update:meas.N;

    pos.noise = nan(3,meas.N);
    % Position update every 10 sample

    pos.sig =  1e-3;   %                     % Sigma position update
    pos.Q = eye(3)*pos.sig^2;
    pos.noise(:,pos.inds_update) = chol(pos.Q)*randn(3, length(pos.inds_update));

    % Release at 15 secs
    pos.noise(:,meas.t > meas.time_point_release) = nan;

    S_true.R = truth.R_nb;
    S_true.w = w_b;
    S_true.p = truth.pos;
    S_true.v = truth.v;
    S_true.omega_dot = truth.v*0; %w_dot_b;

    % A = compute_A_non_center(pos.bodyCorr);  % Appendix A, "Inertial Navigation using an Inertial sensor array"
    % %b_omega_dot_true = -A(1:3,:)*truth.biasAcc;    % adding measurement bias to all omega_dot values
    S_true.b_omega_dot = truth.v*0;%repmat(b_omega_dot_true,1,meas.N);  % omega_dot bias
    b_s_true = -truth.biasAcc;
    S_true.b_s = repmat(b_s_true,1,meas.N);                  % accelerometer bias
    S_true.b_g = repmat(truth.biasGyro,1,meas.N);            % gyro bias is zero by definition as the gyro is absent
    S_true.s = body.acc*0;                              % body linear acceleration in body frame

    settings_default.r = pos.bodyCorr;
    settings_default.T = meas.Ts;
    settings_default.g = settings_default.g * 0;             % in our case g is neglected

    % figure(233);
    % plot(tim, squeeze(body.acc(:,17,:))); 
    % figure(234);
    % plot(tim, body.omega); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END INJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[init, simdata, sensorData, run_settings] = setSimulationParameters(acc, gyro, pos, truth, settings_default);
[resSingle]=runSimulation(init, sensorData, simdata, run_settings, truth, S_true, w_b);
plotResults(meas, resSingle, settings_default.r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getBodyData
function [tim, world, body]=getCFBodyData(dataS, dataParams)
if nargin<2
    dataParams.timfactor = 1;%1e-6;
    dataParams.posfactor = 1;%1e-3;
    dataParams.downsampling = 1;%1000;
end

timfactor = dataParams.timfactor;
if 0
    plot3(dataS.coordx(1,:), dataS.coordy(1,:), dataS.coordz(1,:),'x');axis equal; grid on;
    xlabel('x');ylabel('y');zlabel('z');
    for i=1:length(dataS.coordx(1,:))
        text(dataS.coordx(1,i), dataS.coordy(1,i), dataS.coordz(1,i),num2str(i));
    end
end

indzd=dataS.motor_m1>0;
%  Downsampling
data.coord=[dataS.stateEstimate_x(indzd), dataS.stateEstimate_y(indzd), dataS.stateEstimate_z(indzd)]; 
data.vel  =[dataS.stateEstimate_vx(indzd), dataS.stateEstimate_vy(indzd), dataS.stateEstimate_vz(indzd)];
data.acc  =[dataS.acc_x(indzd),  dataS.acc_y(indzd),  dataS.acc_z(indzd)];
data.gyro =[dataS.gyro_x(indzd), dataS.gyro_y(indzd), dataS.gyro_z(indzd)];
data.angle=[dataS.stateEstimate_yaw, dataS.stateEstimate_pitch, dataS.stateEstimate_roll];

data.tim=dataS.tick(indzd)*timfactor;
tim = data.tim;

world.time=data.tim;
world.pos = data.coord';
world.vel = data.vel';

% calculation of body angular velocity - omega
body.omega = data.gyro';
body.acc = data.acc';
body.yprAngles=data.angle';

disp ''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular_velocitiesMat   -
function [omegaRot] = angular_velocitiesMat(tim, RHist)
skew = @(omega) [0, -omega(3), omega(2);
                 omega(3), 0, -omega(1);
                -omega(2), omega(1), 0];

tlen=length(tim);
dt=zeros(tlen,1);
dt(1:end-1)=diff(tim);  dt(end)=dt(end-1);

dRHist = diff(RHist, 1, 3);  dRHist(:,:,end+1)=dRHist(:,:,end); 

omegaRot = zeros(tlen,3);
for i=1:tlen
    omegaRot_n = squeeze(RHist(:,:,i))'*squeeze(dRHist(:,:,i))/dt(i);
    omegaRot_n = (omegaRot_n-omegaRot_n')/2;
    omegaRot(i,:) = [omegaRot_n(3,2), omegaRot_n(1,3), omegaRot_n(2,1)];
    % 
    % test = expm(skew(omegaRot(i,:))*dt(i))*RHist(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular_velocities    - 2022, Mario García, Angular velocity from Quaternions,
%                          https://mariogc.com/post/angular-velocity-quaternions/
function [omega] = angular_velocities(tim, qHist)
% qHist must be in Body frame
dt=diff(tim);
q1 = qHist(:,1:end-1);
q2 = qHist(:,2:end);
if size(dt,2)==1
    dt=dt';
end
omega = (2 ./ ([1;1;1]*dt)) .* [q1(1,:).*q2(2,:) - q1(2,:).*q2(1,:) - q1(3,:).*q2(4,:) + q1(4,:).*q2(3,:);...
    q1(1,:).*q2(3,:) + q1(2,:).*q2(4,:) - q1(3,:).*q2(1,:) - q1(4,:).*q2(2,:);...
    q1(1,:).*q2(4,:) - q1(2,:).*q2(3,:) + q1(3,:).*q2(2,:) - q1(4,:).*q2(1,:)];

omega(:,end+1) = omega(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular_velocities_w - Angular velocities
% 2002, Schwab, Quaternions, Finite Rotation and Euler Parameters
function [omega] = angular_velocities2(tim, qHist, frame)
% frame: 'w' for world frame, 'b' for body frame
if nargin<3
    frame = 'w';
end
dt=diff(tim);
if size(dt,2)==1
    dt=dt';
end
qdot = diff(qHist,1,2)./([1,1,1,1]'*dt);  qdot(:,end+1)=qdot(:,end);
omega = zeros(4,length(qdot(1,:)));

if frame=='w'
    rotmat = QmatbarT(qHist);
else
    rotmat = QmatT(qHist);
end

for i=1:length(tim)
    omega(:,i) = 2*rotmat(:,:,i)*qdot(:,i);  %2*QmatT(qHist)*qdot';
end
omega = ([1;1;1]*sign(omega(1,:))).*omega(2:4,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular_accelerations   - 2002, Schwab, Quaternions, Finite Rotation and Euler Parameters
% 2002, Schwab, Quaternions, Finite Rotation and Euler Parameters
function [omegadot] = angular_accelerations(tim, qHist, frame)
% frame: 'w' for world frame, 'b' for body frame
if nargin<3
    frame = 'w';
end
dt=diff(tim);
qdot = diff(qHist,1,2)./([1;1;1;1]*dt');    qdot(:,end+1)=qdot(:,end);
qdotdot = diff(qdot,1,2)./([1;1;1;1]*dt');  qdotdot(:,end+1)=qdotdot(:,end);

if frame=='w'
    rotmat = QmatbarT(qHist);
else
    rotmat = QmatT(qHist);
end

omegadot = zeros(4, length(tim));
for i=1:length(tim)   %    2*QmatT(qHist)*qdotdot'+2*[norm(qdot)^2; 0; 0; 0];
    omegadot(:,i) = 2*rotmat(:,:,i)*qdotdot(:,i)+2*[norm(qdot(:,i))^2; 0; 0; 0];
end
omegadot = ([1;1;1]*sign(omegadot(1,:))).*omegadot(2:4,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Qmat
function [qmat] = Qmat(q)
qmat = zeros(4,4,length(q(1,:)));
qmat(1,1,:) = q(1,:);  qmat(1,2,:) = -q(2,:); qmat(1,3,:) = -q(3,:); qmat(1,4,:) = -q(4,:);
qmat(2,1,:) = q(2,:);  qmat(2,2,:) =  q(1,:); qmat(2,3,:) = -q(4,:); qmat(2,4,:) =  q(3,:);
qmat(3,1,:) = q(3,:);  qmat(3,2,:) =  q(4,:); qmat(3,3,:) =  q(1,:); qmat(3,4,:) = -q(2,:);
qmat(4,1,:) = q(4,:);  qmat(4,2,:) = -q(3,:); qmat(4,3,:) =  q(2,:); qmat(4,4,:) =  q(1,:);

% qmat = [q(1), -q(2), -q(3), -q(4);
%         q(2),  q(1), -q(4),  q(3);
%         q(3),  q(4),  q(1), -q(2);
%         q(4), -q(3),  q(2),  q(1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QmatT
function [qmatT] = QmatT(q)
qmatT = zeros(4,4,length(q(1,:)));
qmatT(1,1,:) =  q(1,:);  qmatT(1,2,:) =  q(2,:); qmatT(1,3,:) =  q(3,:); qmatT(1,4,:) =  q(4,:);
qmatT(2,1,:) = -q(2,:);  qmatT(2,2,:) =  q(1,:); qmatT(2,3,:) =  q(4,:); qmatT(2,4,:) = -q(3,:);
qmatT(3,1,:) = -q(3,:);  qmatT(3,2,:) = -q(4,:); qmatT(3,3,:) =  q(1,:); qmatT(3,4,:) =  q(2,:);
qmatT(4,1,:) = -q(4,:);  qmatT(4,2,:) =  q(3,:); qmatT(4,3,:) = -q(2,:); qmatT(4,4,:) =  q(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Qmatbar
function [qmatbar] = Qmatbar(q)
qmatbar = zeros(4,4,length(q(1,:)));
qmatbar(1,1,:) = q(1,:);  qmatbar(1,2,:) = -q(2,:); qmatbar(1,3,:) = -q(3,:); qmatbar(1,4,:) = -q(4,:);
qmatbar(2,1,:) = q(2,:);  qmatbar(2,2,:) =  q(1,:); qmatbar(2,3,:) =  q(4,:); qmatbar(2,4,:) = -q(3,:);
qmatbar(3,1,:) = q(3,:);  qmatbar(3,2,:) = -q(4,:); qmatbar(3,3,:) =  q(1,:); qmatbar(3,4,:) =  q(2,:);
qmatbar(4,1,:) = q(4,:);  qmatbar(4,2,:) =  q(3,:); qmatbar(4,3,:) = -q(2,:); qmatbar(4,4,:) =  q(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QmatbarT
function [qmatbarT] = QmatbarT(q)
qmatbarT = zeros(4,4,length(q(1,:)));
qmatbarT(1,1,:) =  q(1,:);  qmatbarT(1,2,:) =  q(2,:); qmatbarT(1,3,:) =  q(3,:); qmatbarT(1,4,:) =  q(4,:);
qmatbarT(2,1,:) = -q(2,:);  qmatbarT(2,2,:) =  q(1,:); qmatbarT(2,3,:) = -q(4,:); qmatbarT(2,4,:) =  q(3,:);
qmatbarT(3,1,:) = -q(3,:);  qmatbarT(3,2,:) =  q(4,:); qmatbarT(3,3,:) =  q(1,:); qmatbarT(3,4,:) = -q(2,:);
qmatbarT(4,1,:) = -q(4,:);  qmatbarT(4,2,:) = -q(3,:); qmatbarT(4,3,:) =  q(2,:); qmatbarT(4,4,:) =  q(1,:);
% qmatbar = [q(1), -q(2), -q(3), -q(4);
%            q(2),  q(1),  q(4), -q(3);
%            q(3), -q(4),  q(1),  q(2);
%            q(4),  q(3), -q(2),  q(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initSamplingData
function [meas, acc, gyro, pos]=initSamplingData()
meas.Fs = 1000;                  % [Hz] Sampling frequency
meas.Ts = 1/meas.Fs;              % [s] Sampling period
meas.time_factor=1;
meas.T_end = 1;                 % [s] End time of simulations
meas.t = 0:meas.Ts:meas.T_end - meas.Ts;
meas.N = length(meas.t);
meas.noiseFactor=1; %1e-4;
meas.biasFactor=1;

F_pos = 1000;                     % [Hz] Time frequency of position update in simulation
assert(mod(meas.Fs,F_pos) == 0)
pos.N_update = meas.Fs/F_pos;     % update position every N steps in simulation
pos.inds_update = 1:pos.N_update:meas.N;

% Low dynamics
% N_periods = 1.5 % around 2000 deg/s in norm(w)
% High dynamics

meas.time_point_release = 40;

gyroNoiseSigContDeg = meas.noiseFactor*1/sqrt(500);   % Gyro continous noise
accNoiseSigCont = meas.noiseFactor*100/sqrt(500)*(1e-8);     % Acc continous noise

gyroNoiseSigCont = deg2rad(gyroNoiseSigContDeg);

acc.noiseSigDisc = sqrt(meas.Fs)*accNoiseSigCont;     % Discrete Acc noise
gyro.noiseSigDisc = sqrt(meas.Fs)*gyroNoiseSigCont;   % Discrete gyro noise
rad2deg(gyro.noiseSigDisc);

pos.sig = 1e-1;        % Sigma position update

% estimate a constant bias
acc.Q_b_a = meas.biasFactor*0;                      % accelerometer bias value, From  Carlsson2021
gyro.Q_b_g = meas.biasFactor*0;                     % gyro bias value, From  Carlsson2021
acc.init_b_a = meas.noiseFactor*0.2*1e-8;               % initial acceleration bias
gyro.init_b_g = meas.noiseFactor*deg2rad(1)*0;       % initial gyro bias

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% genBodyFrameAngleTrajectory
function [phi, theta, vels]=genBodyFrameAngleTrajectory(meas)
global showfigures;

% phi in [0, pi]
inp.phi = "sinus";
params = struct;
params.A = 5;
params.b = 1;
params.f = 10/6.;% Three periods
inp.phi_params = params;

% theta in [0, 2pi]
inp.theta = "sinus";
params = struct;
params.A = 10;
params.b = 1;
params.f = 20/6.; % two periods
inp.theta_params = params;

m = get_spherical_motion(meas.t, inp);   % calculate phi,theta, phi_dot, theta_dot history

vels = compute_omega(m); % calculation of omega and omega_dot in navigation and body frames

% f = figure(1);
% set(f,'WindowStyle','docked')
% clf

if showfigures
    figure(1);
    vels = mergeStruct(m, vels);
    plot_angular_velocities(meas.t, vels)
    scale_figure(gcf, 1.5);
    figure(2)
    hold on
    plot(meas.t, rad2deg(norm_time(vels.w_b)), "b-")
    plot(meas.t, rad2deg(norm_time(vels.w_n)), "r--")
end

phi = m.phi;
theta = m.theta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getNav2BodyRotHistory
function [truth, w_b, w_dot_b] = getNav2BodyRotHistory(phi, theta, meas, vels)
global showfigures;

truth.R_bn = R_bn_spherical(phi,theta, vels.Rb_om, vels.Rn_om); % calculate navigation to body frame rotation matrix history
R_offset = eye(3);

w_b = R_offset*vels.w_b;
w_dot_b = R_offset*vels.w_dot_b;


% Make initial rotation eye(3)
for n = 1:meas.N
    truth.R_bn(:,:,n) = R_offset'*truth.R_bn(:,:,n);
end
truth.R_nb = zeros(size(truth.R_bn));
for n = 1:meas.N
    truth.R_nb(:,:,n) = truth.R_bn(:,:,n)';
end
%%

if showfigures
    figure(3);
    eul = my_rotm2eul(truth.R_nb);
    plot(meas.t,rad2deg(eul'))
    legend("roll", "pitch", "yaw")
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getSensorPositions
function [pos]=getSensorPositions(pos)
global showfigures;

% sensor positions in body frame
pos.body = [
    0    0    0
    0    0.1  0
    1    0    0
    0    0   0.1]';
% i_IMU_select = [3,5,27,29];
% i_IMU_select = sort([i_IMU_select, i_IMU_select+1])
i_IMU_select = 1:4; % select specific IMUs from the list in case that not all IMU are used
pos.body = pos.body(:,i_IMU_select);   % Select specific IMUs for calculations
%%
pos.bodyCorr = pos.body - mean(pos.body,2);   % corrected by centroid of sensor positions

%%
if showfigures
    figure(4);
    subplot(1,2,1)
    plot(pos.bodyCorr(1,1:2:end)*1e3, pos.bodyCorr(2,1:2:end)*1e3, "x")
    axis("equal")
    grid on
    L = 15;
    % xlim([-L,L])
    % ylim([-L,L])
    % xlim([-5,5]*1e-3)
    xlabel("[mm]")
    ylabel("[mm]")
    title("Topside")
    subplot(1,2,2)
    plot(pos.bodyCorr(1,2:2:end)*1e3, pos.bodyCorr(2,2:2:end)*1e3, "x")
    grid on
    axis("equal")
    % xlim([-L,L])
    % ylim([-L,L])
    xlabel("[mm]")
    ylabel("[mm]")
    title("underside")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getGeneralSettings
function [settings_default, inds]=getGeneralSettings(meas, pos)
K = size(pos.body,2); % Number of acc triads
inds = reshape(1:3*K,3,[]);

settings_default = struct;
settings_default.save_full_covariances = false;
settings_default.save_jacobians = false;
settings_default.verbose = true;
settings_default.T = meas.Ts;                            % Sampling period
settings_default.g = [0; 0; -9.81];

settings_default.r = pos.bodyCorr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generateMotion
function [truth] = generateMotionTrajectory(meas, truth)
global showfigures;

f_pos = [13 1 150]/meas.T_end*meas.time_factor; % frequency for each axis


truth.pos = zeros(3,meas.N);
truth.v = zeros(3,meas.N);
truth.a = zeros(3,meas.N);

for i = 1:3
    f_i = f_pos(i);
    truth.pos(i,:) = sin(f_i*meas.t);
    truth.v(i,:) = f_i*cos(f_i*meas.t);
    truth.a(i,:) = -f_i^2*sin(f_i*meas.t);
end
%%

if showfigures
    figure(5);
    subplot(2,1,1)
    hold on
    for i = 1:3
        plot(meas.t, truth.pos(i,:))

    end

    grid on
    legend("x","y","z")
    ylabel("[m]")
    xlabel("Time [s]")
    title("True Navigation Position");
    subplot(2,1,2)
    plot3(truth.pos(1,:),truth.pos(2,:),truth.pos(3,:))
    grid on
    xlabel("x[m]")
    ylabel("y[m]")
    zlabel("z[m]")
    axis equal


    figure(6)
    hold on
    for i = 1:3
        plot(meas.t, truth.v(i,:))

    end
    grid on
    legend("x","y","z")
    ylabel("[m/s]")
    xlabel("Time [s]")
    title("True Navigation Velocity");


    figure(7)
    hold on
    for i = 1:3
        plot(meas.t, truth.a(i,:))

    end
    grid on
    legend("x","y","z")
    ylabel("[m/s]")
    xlabel("Time [s]")
    title("True Navigation Acceleration");

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% addNoiseAndBias
function [S_true, acc, gyro, pos, truth] = addNoiseAndBias(meas, pos, settings_default, gyro, acc, w_b, w_dot_b, truth, inds)
% Zero acceleration, stay at the same position. Measure gravity in body frame

K = size(pos.body,2); % Number of acc triads
truth.acc_u = zeros(3*K,meas.N);
rot_comp = zeros(3*K,meas.N);
trans_comp = zeros(3,meas.N);
for n = 1:meas.N
    W = skew_sym(w_b(:,n))^2 + skew_sym(w_dot_b(:,n));
    trans_comp(:,n) = truth.R_bn(:,:,n)*(truth.a(:,n) - settings_default.g);   % convert body center accelerations
    % from inertial frame to body frame
    for k = 1:K
        kk = inds(:,k);
        r_k = pos.bodyCorr(:,k);                               % sensors position
        rot_comp(kk,n) = W*r_k;                                % rotational acceleration component
        truth.acc_u(kk,n) = trans_comp(:,n) + rot_comp(kk,n);  % total acceleration in body frame
    end
end

truth.gyro_u = w_b;                                            % Angular velocity in body coordinates
%%
rng("default");  % set random number generator type

truth.biasGyro = gyro.init_b_g*randn(3,1)/sqrt(K);
truth.biasAcc = acc.init_b_a*randn(3*K,1);

A = compute_A_non_center(pos.bodyCorr);  % Appendix A, "Inertial Navigation using an Inertial sensor array"
b_omega_dot_true = -A(1:3,:)*truth.biasAcc;    % adding measurement bias to all omega_dot values
b_s_true = -A(4:6,:)*truth.biasAcc;

S_true = struct;
S_true.R = truth.R_nb;
S_true.w = w_b;
S_true.p = truth.pos;
S_true.v = truth.v;
S_true.omega_dot = w_dot_b;

S_true.b_omega_dot = repmat(b_omega_dot_true,1,meas.N);  % omega_dot bias
S_true.b_s = repmat(b_s_true,1,meas.N);                  % accelerometer bias
S_true.b_g = repmat(truth.biasGyro,1,meas.N);            % gyro bias
S_true.s = trans_comp;                              % body linear acceleration in body frame

acc.Q = eye(3*K)*acc.noiseSigDisc^2;
gyro.Q = eye(3)*gyro.noiseSigDisc^2/K;
pos.Q = eye(3)*pos.sig^2;

acc.noise = chol(acc.Q)*randn(3*K,meas.N);
gyro.noise = chol(gyro.Q)*randn(3,meas.N);
pos.noise = nan(3,meas.N);
% Position update every 10 sample
pos.noise(:,pos.inds_update) = chol(pos.Q)*randn(3, length(pos.inds_update));

% Release at 15 secs
pos.noise(:,meas.t > meas.time_point_release) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setSimulationParameters
function [init, simdata, sensorData, run_settings] = setSimulationParameters(acc, gyro, pos, truth, settings_default)
% Initial values

rng("default");
K = 1; % Number of acc triads

init = struct;
init.cov.R = eye(3)*deg2rad(1)^2;
init.cov.p = pos.Q*2;
init.cov.v = eye(3)*0.1^2;
init.cov.omega = gyro.Q*2;

init.mean.R = truth.R_nb(:,:,1)*expSO3(chol(init.cov.R)*randn(3,1));
init.mean.p = truth.pos(:,1) ;
init.mean.v = truth.v(:,1) + chol(init.cov.v)*randn(3,1);
init.mean.omega = truth.gyro_u(:,1);

init.mean.b_a = zeros(3*K,1);
init.cov.b_a = acc.init_b_a^2*eye(3*K)*3;

init.mean.b_g = zeros(3,1);
init.cov.b_g = eye(3)*gyro.init_b_g^2/K*3;

% Run Filters
 
% error("asd")
% Models

simdata = struct;
simdata.gyroscope_1st_order = settings_default;
simdata.gyroscope_1st_order.input_accelerometers = true;
simdata.gyroscope_1st_order.propagate_bias_s = true;
simdata.gyroscope_1st_order.propagate_bias_gyro = true;
simdata.gyroscope_1st_order.propagate_position = true;
simdata.gyroscope_1st_order.propagate_velocity = true;
simdata.gyroscope_1st_order.get_model = @D_LG_EKF_Gyro_1st_v4;
simdata.gyroscope_1st_order.r = pos.bodyCorr;
simdata.gyroscope_1st_order.N_a = K;
simdata.gyroscope_1st_order.Q_acc = acc.noiseSigDisc^2*eye(3*K);
simdata.gyroscope_1st_order.Q_bias_acc = acc.Q_b_a^2*eye(3*K);
simdata.gyroscope_1st_order.R_pos = pos.sig^2*eye(3);
simdata.gyroscope_1st_order.Q_gyro = gyro.noiseSigDisc^2*eye(3)/K;
simdata.gyroscope_1st_order.Q_bias_gyro = gyro.Q_b_g^2*eye(3)/K;
simdata.gyroscope_1st_order.do_gyro_updates = false;
simdata.gyroscope_1st_order.do_position_updates = true;
simdata.gyroscope_1st_order.do_rotation_updates = false;
simdata.gyroscope_1st_order.label = "1st order gyroscope";

sensorData = struct;
sensorData.acc_measurements = truth.acc_u;
sensorData.gyro_measurements = truth.gyro_u;
sensorData.position_measurements = truth.pos;

% Run filters
run_settings = struct;
run_settings.compute_error = true;
run_settings.verbose = true;
run_settings.save_input = true;
run_settings.save_residuals = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% runSimulation
function [resSingle]=runSimulation(init, sensorData, simdata, run_settings, truth, S_true, w_b)
resSingle = struct;
resSingle.label = "";
[~,resSingle.gyroscope_1st_order] = run_filter(sensorData, init, simdata.gyroscope_1st_order,run_settings, S_true);

%
true_traj = struct;
true_traj.filt.mean.b_g = S_true.b_g;
true_traj.filt.mean.b_omega_dot = S_true.b_omega_dot;
true_traj.filt.mean.b_s = S_true.b_s;
true_traj.filt.mean.p = truth.pos;
true_traj.filt.mean.v = truth.v;
true_traj.filt.mean.R = truth.R_nb;
true_traj.filt.mean.w = w_b;
true_traj.filt.std = struct;
true_traj.label = "True";

resSingle.True = struct;
resSingle.True.res = true_traj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotcoord
function plotcoord(p,R, linewidth) 
quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), 'r','linewidth',linewidth); hold on; axis equal; grid on; 
quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), 'g','linewidth',linewidth), 
quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), 'b','linewidth',linewidth);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotResults
function plotResults(meas, resSingle, r)

data1=resSingle.gyroscope_1st_order.res.filt;
dataT=resSingle.True.res.filt;

% plot(meas.t, dataT.mean.p); hold on; plot(meas.t, data1.mean.p); 

ntim=length(meas.t);
for i=1:ntim
    quat1(:,i) = dcm2quat(squeeze(data1.mean.R(:,:,i)));
    quatT(:,i) = dcm2quat(squeeze(dataT.mean.R(:,:,i)));
    if i>1
        quat1(:,i)=quat1(:,i)*sign(dot(quat1(:,i),quat1(:,i-1)));
        quatT(:,i)=quatT(:,i)*sign(dot(quatT(:,i),quatT(:,i-1)));
    end
end
% plot(meas.t, quatT); hold on; plot(meas.t, quat1); 
for i=1:round(ntim/100):ntim
    p1 = data1.mean.p(:,i);  R1 = data1.mean.R(:,:,i);
    pT = dataT.mean.p(:,i);  RT = dataT.mean.R(:,:,i);
    plotcoord (pT, RT,1);
    plotcoord (p2, R2,3);
end

plotBias = 0;
plotRotation = 1;
plotPosition = 1;

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
if plotBias
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order";
        "True"];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "filt"];
    [fig,~] = plot_bias_s(resSingle, meas.t, cases_plot, plotOpt_k, "normal", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"];
    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    [fig,~] = plot_bias_s(resSingle, meas.t, cases_plot, plotOpt_k, "error", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    extras.show_mean = false;
    [fig,~] = plot_bias_s(resSingle, meas.t, cases_plot, plotOpt_k, "error_log", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order"; "True"];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "filt"];
    [fig,~] = plot_bias_omega_dot(resSingle, meas.t, cases_plot, plotOpt_k, "normal", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    extras.show_mean = false;
    [fig,~] = plot_bias_omega_dot(resSingle, meas.t, cases_plot, plotOpt_k, "error_log", extras);
    scale_figure(fig,1.5);

    % Bias gyroscope
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"; "True"];
    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "filt"];
    [fig,~] = plot_bias_gyroscopes(resSingle, meas.t, cases_plot, plotOpt_k, "normal", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"];
    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    extras.show_mean = false;
    [fig,~] = plot_bias_gyroscopes(resSingle, meas.t, cases_plot, plotOpt_k, "error_log", extras);
    scale_figure(fig,1.5);
end
%
if plotRotation
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"];
    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    [fig,~] = plot_rotation(resSingle, meas.t, cases_plot, plotOpt_k, "error", extras);
    scale_figure(fig,1.5);
end
%
if plotPosition
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order"; "True"];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "filt"];
    [fig,~] = plot_navigation_position(resSingle, meas.t, cases_plot, plotOpt_k, "normal", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order";];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    [fig,~] = plot_navigation_position(resSingle, meas.t, cases_plot, plotOpt_k, "error", extras);
    scale_figure(fig,1.5);
    %
    cases_plot = [ "accelerometer_array_2nd_order";
        "accelerometer_array_1st_order";
        "gyroscope_2nd_order";
        "gyroscope_1st_order";];

    plotOpt_k = struct;
    extras = struct;
    extras.path = ["res", "err"];
    extras.show_mean = false;
    [fig,~] = plot_navigation_position(resSingle, meas.t, cases_plot, plotOpt_k, "error_log", extras);
    scale_figure(fig,1.5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% run_filter
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% compute_error_tot
function err = compute_error_tot(S_tot, S_ref)

err = struct;
err.mean = compute_error(S_tot.filt.mean, S_ref);

% Copy standard deviations
fields = fieldnames(err.mean);
for idx = 1:length(fields)
    err.std.(fields{idx}) = S_tot.filt.std.(fields{idx});
end



