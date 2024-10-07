%% Single Trajectory simulation 
%% 
function run_script_single_3_LSdyna
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
    data=load('../../Data/lssimDatad10.mat');
    [tim, world, body, npoints]=getBodyData(data);
    truth.pos = world.origin.pos;
    truth.v = world.origin.vel;
    truth.a = world.origin.acc;
    truth.R_nb = world.rot_wb;
    truth.R_bn = world.rot_bw;

    ntim=length(tim);
    acc_u=zeros(3*npoints,ntim);
    body.acc_origin = body.origin.acc;
    for i=1:ntim
        bacc=body.acc(:,:,i);
        acc_u(:,i) = bacc(:);
    end
    truth.acc_u=acc_u;
    truth.gyro_u=truth.v*0;

    meas.N=ntim;
    meas.t=tim;
    meas.T_end = tim(end);
    meas.Ts=tim(2)-tim(1);
    meas.Fs=1/meas.Ts;

    truth.biasAcc = acc.init_b_a*randn(3*npoints,1);
    acc.Q = eye(3*npoints)*acc.noiseSigDisc^2;
    acc.noise = chol(acc.Q)*randn(3*npoints,meas.N);
    gyro.noise = chol(gyro.Q)*randn(3,meas.N);

    F_pos = meas.Fs/1;                     % [Hz] Time frequency of position update in simulation, must be divider of sampling frequency
    % assert(mod(meas.Fs,F_pos) == 0)
    pos.N_update = floor(meas.Fs/F_pos);     % update position every N steps in simulation
    pos.inds_update = 1:pos.N_update:meas.N;

    pos.body = [data.coordx(1,:); data.coordx(2,:); data.coordx(3,:)];
    pos.bodyCorr = squeeze(body.pos(:,:,1));
    pos.noise = nan(3,meas.N);
    % Position update every 10 sample
    pos.noise(:,pos.inds_update) = chol(pos.Q)*randn(3, length(pos.inds_update));

    % Release at 15 secs
    pos.noise(:,meas.t > meas.time_point_release) = nan;

    S_true.R = truth.R_nb;
    S_true.w = w_b;
    S_true.p = truth.pos;
    S_true.v = truth.v;
    S_true.omega_dot = truth.v*0;%w_dot_b;

    A = compute_A_non_center(pos.bodyCorr);  % Appendix A, "Inertial Navigation using an Inertial sensor array"
    %b_omega_dot_true = -A(1:3,:)*truth.biasAcc;    % adding measurement bias to all omega_dot values
    S_true.b_omega_dot = truth.v*0;%repmat(b_omega_dot_true,1,meas.N);  % omega_dot bias 
    b_s_true = -A(4:6,:)*truth.biasAcc;
    S_true.b_s = repmat(b_s_true,1,meas.N);                  % accelerometer bias
    S_true.b_g = repmat(truth.biasGyro,1,meas.N);            % gyro bias 
    S_true.s = body.acc_origin;                              % body linear acceleration in body frame

    settings_default.r = pos.bodyCorr;
    settings_default.T = meas.Ts;

    figure(233);
    plot(tim, squeeze(body.acc(1,1,:))); hold on; plot(tim, squeeze(body.acc(2,1,:))); plot(tim, squeeze(body.acc(3,1,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END INJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[init, simdata, sensorData, run_settings] = setSimulationParameters(acc, gyro, pos, truth, settings_default);
[resSingle]=runSimulation(init, sensorData, simdata, run_settings, truth, S_true, w_b);
plotResults(meas, resSingle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getBodyData
function [tim, world, body, npoints]=getBodyData(data)
% cg.x=median(data.coordx(:,:), 2);
% cg.y=cg.x*0;
% cg.z=median(data.coordz(:,:), 2);
indz=data.tim==0;
data.disx(indz,:)=[];
data.disy(indz,:)=[];
data.disz(indz,:)=[];
data.velx(indz,:)=[];
data.vely(indz,:)=[];
data.velz(indz,:)=[];
data.accx(indz,:)=[];
data.accy(indz,:)=[];
data.accz(indz,:)=[];
data.coordx(indz,:)=[];
data.coordy(indz,:)=[];
data.coordz(indz,:)=[];
data.tim(indz)=[];
downsampling = 100;
indzd = 1:downsampling:length(data.tim);

timfactor = 1e-3;
posfactor = 1e-3;
velfactor = posfactor/timfactor;
accfactor = posfactor/timfactor/timfactor;

data.disx=data.disx(indzd,:)*posfactor; data.disy=data.disy(indzd,:)*posfactor; data.disz=data.disz(indzd,:)*posfactor;
data.velx=data.velx(indzd,:)*velfactor; data.vely=data.vely(indzd,:)*velfactor; data.velz=data.velz(indzd,:)*velfactor;
data.accx=data.accx(indzd,:)*accfactor; data.accy=data.accy(indzd,:)*accfactor; data.accz=data.accz(indzd,:)*accfactor;
data.coordx=data.coordx(indzd,:)*posfactor; data.coordy=data.coordy(indzd,:)*posfactor; data.coordz=data.coordz(indzd,:)*posfactor;
data.tim=data.tim(indzd)*timfactor;
tim = data.tim;

dt=diff(tim); dt(end+1)=dt(end);

if downsampling>1
    data.velx = diff(data.coordx)./diff(data.tim);    data.vely = diff(data.coordy)./diff(data.tim);   data.velz = diff(data.coordz)./diff(data.tim);
    data.velx = [data.velx; data.velx(end,:)];      data.vely = [data.vely; data.vely(end,:)];       data.velz = [data.velz; data.velz(end,:)];

    data.accx = diff(data.velx)./diff(data.tim);    data.accy = diff(data.vely)./diff(data.tim);   data.accz = diff(data.velz)./diff(data.tim);
    data.accx = [data.accx; data.accx(end,:)];      data.accy = [data.accy; data.accy(end,:)];       data.accz = [data.accz; data.accz(end,:)];

    data.disx = cumsum(data.velx.*dt,1);  data.disy = cumsum(data.vely.*dt,1);   data.disz = cumsum(data.velz.*dt,1);
end

npoints = length(data.coordx(1,:));
ntim = length(data.tim);

dirX=[data.coordx(:,15)-data.coordx(:,16), data.coordy(:,15)-data.coordy(:,16), data.coordz(:,15)-data.coordz(:,16)];
d24_16=[data.coordx(:,24)-data.coordx(:,16), data.coordy(:,24)-data.coordy(:,16), data.coordz(:,24)-data.coordz(:,16)];
for i=1:length(tim)
    dirZ(i,:)=cross(dirX(i,:), cross(d24_16(i,:), dirX(i,:)));
    dirY(i,:)=cross(dirX(i,:), dirZ(i,:));
    dirX(i,:)=dirX(i,:)/norm(dirX(i,:));
    dirY(i,:)=dirY(i,:)/norm(dirY(i,:));
    dirZ(i,:)=dirZ(i,:)/norm(dirZ(i,:));
    disp ''
end
world.frame.x=dirX;
world.frame.y=dirY;
world.frame.z=dirZ;

world.origin.pos = [data.coordx(:,15), data.coordy(:,15), data.coordz(:,15)]';
world.origin.vel = [data.velx(:,15), data.vely(:,15), data.velz(:,15)]';
world.origin.acc = [data.accx(:,15), data.accy(:,15), data.accz(:,15)]';

world.rot_wb = zeros(3,3,ntim);
world.rot_bw = zeros(3,3,ntim);
for i=1:length(tim)
    world.rot_wb(:,:,i) = [dirX(i,:)', dirY(i,:)', dirZ(i,:)'];
    world.rot_bw(:,:,i) = [dirX(i,:)', dirY(i,:)', dirZ(i,:)']';
end
world.omegaMat = angular_velocitiesMat(data.tim, world.rot_wb);
% world.omega = angular_velocities(data.tim, world.quat_wb);
% world.omega_2 = angular_velocities2(data.tim, world.quat_wb);
% world.omegadot = angular_accelerations(data.tim, world.quat_wb);

body.pos = zeros(3,npoints, length(tim));         body.vel = zeros(3,npoints, length(tim));         body.acc = zeros(3,npoints, length(tim)); 
body.origin.pos = zeros(3, length(tim));  body.origin.vel = zeros(3, length(tim));  body.origin.acc = zeros(3, length(tim));
world.pos = zeros(3,npoints, length(tim)); 
for i=1:length(tim)
    rot_bw = squeeze(world.rot_wb(:,:,i));
    world.pos(:,:,i)=[data.coordx(i,:)-world.origin.pos(1,1); 
                      data.coordy(i,:)-world.origin.pos(1,2);
                      data.coordz(i,:)-world.origin.pos(1,3)];
    world.vel(:,:,i)=[data.velx(i,:); data.vely(i,:); data.velz(i,:)];
    world.acc(:,:,i)=[data.accx(i,:); data.accy(i,:); data.accz(i,:)];

    body.pos(:,:,i) = rot_bw*world.pos(:,:,i);              body.vel(:,:,i) = rot_bw*world.vel(:,:,i);              body.acc(:,:,i) = rot_bw*world.acc(:,:,i);
    body.origin.pos(:,i) = rot_bw*world.origin.pos(:,i);    body.origin.vel(:,i) = rot_bw*world.origin.vel(:,i);    body.origin.acc(:,i) = rot_bw*world.origin.acc(:,i);
end
dpos = diff(world.pos,1,3);  dpos(:,:,end+1)=dpos(:,:,end);
ddpos = diff(dpos,1,3);  ddpos(:,:,end+1)=ddpos(:,:,end);
for i=1:length(tim)
    world.vel2(:,:,i)=dpos(:,:,i)./dt(i);
    world.acc2(:,:,i)=ddpos(:,:,i)./dt(i)/dt(i);
end
% plot3(squeeze(world.pos(1,:,end)), squeeze(world.pos(2,:,end)), squeeze(world.pos(3,:,end)), 'x'); grid on; axis equal; hold on
% plot3(squeeze(world.pos(1,:,1)), squeeze(world.pos(2,:,1)), squeeze(world.pos(3,:,1)), 'x'); grid on; axis equal

% figure(354)
% plot(tim, squeeze(world.vel(1,1,:)), tim, squeeze(world.vel(2,1,:)), tim, squeeze(world.vel(3,1,:))); grid on; hold on
% % plot(tim, squeeze(world.vel2(1,1,:)), tim, squeeze(world.vel2(2,1,:)), tim, squeeze(world.vel2(3,1,:))); grid on; hold on
% plot(tim, sqrt(squeeze(world.vel(1,1,:)).^2+squeeze(world.vel(2,1,:)).^2+squeeze(world.vel(3,1,:)).^2));
% % plot(tim, sqrt(squeeze(world.vel2(1,1,:)).^2+squeeze(world.vel2(2,1,:)).^2+squeeze(world.vel2(3,1,:)).^2));
% legend('x', 'y', 'z', 'x2', 'y2', 'z2', 'norm1', 'norm2')
% figure(435)
% plot(tim, squeeze(world.acc(1,1,:)), tim, squeeze(world.acc(2,1,:)), tim, squeeze(world.acc(3,1,:))); grid on; hold on
% plot(tim, squeeze(world.acc2(1,1,:)), tim, squeeze(world.acc2(2,1,:)), tim, squeeze(world.acc2(3,1,:))); grid on; hold on
% plot(tim, sqrt(squeeze(world.acc(1,1,:)).^2+squeeze(world.acc(2,1,:)).^2+squeeze(world.acc(3,1,:)).^2));
% plot(tim, sqrt(squeeze(world.acc2(1,1,:)).^2+squeeze(world.acc2(2,1,:)).^2+squeeze(world.acc2(3,1,:)).^2));
% legend('x', 'y', 'z', 'x2', 'y2', 'z2', 'norm1', 'norm2')

% figure(87);
% subplot(2,1,1);
% plot(tim, squeeze(body.acc(1,1,:))); hold on; plot(tim, squeeze(body.acc(2,1,:))); plot(tim, squeeze(body.acc(3,1,:)));
% legend('x', 'y', 'z');
% xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
% title('Body Acceleration');
% subplot(2,1,2);
% plot(tim, squeeze(world.acc(1,1,:))); hold on; plot(tim, squeeze(world.acc(2,1,:))); plot(tim, squeeze(world.acc(3,1,:)));
% legend('x', 'y', 'z');
% xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
% title('World Acceleration');

% plot(tim, squeeze(body.vel(1,1,:))); hold on; plot(tim, squeeze(body.vel(2,1,:))); plot(tim, squeeze(body.vel(3,1,:)));
% legend('x', 'y', 'z');
% xlabel('Time (s)'); ylabel('Velocity (m/s)');
% title('Body Velocity');
% subplot(3,1,3);
% figure(45);
% subplot(2,1,1)
% plot(tim, world.origin.pos(:,1)); hold on; grid on;plot(tim, world.origin.pos(:,2)); plot(tim, world.origin.pos(:,3)); 
% legend('x', 'y', 'z');
% xlabel('Time (s)'); ylabel('Position (m)');
% title('World Position');
% subplot(2,1,2)
% plot(tim, squeeze(body.pos(1,1,:))); hold on; grid on; plot(tim, squeeze(body.pos(2,1,:))); plot(tim, squeeze(body.pos(3,1,:)));
% legend('x', 'y', 'z');
% xlabel('Time (s)'); ylabel('Position (m)');
% title('Body Position');

disp ''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% angular_velocitiesMat   - 
function [omegaRot] = angular_velocitiesMat(tim, RHist)

tlen=length(tim);
dt=zeros(tlen,1);
dt(1:end-1)=diff(tim);  dt(end)=dt(end-1);
dRHist = zeros(3,3,tlen);
dRHist(:,:,1:end-1) = diff(RHist, 1, 3);

omegaRot = zeros(tlen,3);
for i=1:tlen
    omegaRot_n = squeeze(RHist(:,:,i))'*squeeze(dRHist(:,:,i))/dt(i);
    omegaRot_n = (omegaRot_n-omegaRot_n')/2;
    omegaRot(i,:) = [omegaRot_n(2,3), omegaRot_n(3,1), omegaRot_n(1,2)];
end

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
accNoiseSigCont = meas.noiseFactor*100/sqrt(500);     % Acc continous noise 

gyroNoiseSigCont = deg2rad(gyroNoiseSigContDeg); 

acc.noiseSigDisc = sqrt(meas.Fs)*accNoiseSigCont;     % Discrete Acc noise
gyro.noiseSigDisc = sqrt(meas.Fs)*gyroNoiseSigCont;   % Discrete gyro noise
rad2deg(gyro.noiseSigDisc);

pos.sig = 1e-1;        % Sigma position update

% estimate a constant bias
acc.Q_b_a = meas.biasFactor*0;                      % accelerometer bias value, From  Carlsson2021
gyro.Q_b_g = meas.biasFactor*0;                     % gyro bias value, From  Carlsson2021
acc.init_b_a = meas.noiseFactor*0.2;               % initial acceleration bias
gyro.init_b_g = meas.noiseFactor*deg2rad(1);       % initial gyro bias

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
K = size(pos.body,2); % Number of acc triads

init = struct;
init.cov.R = eye(3)*deg2rad(1)^2;
init.cov.p = pos.Q*2; 
init.cov.v = eye(3)*0.1^2;
init.cov.omega = gyro.Q*2;

init.mean.R = truth.R_nb(:,:,1)*expSO3(chol(init.cov.R)*randn(3,1));
init.mean.p = truth.pos(:,1) + pos.noise(:,1); 
init.mean.v = truth.v(:,1) + chol(init.cov.v)*randn(3,1);
init.mean.omega = truth.gyro_u(:,1) + gyro.noise(:,1);

init.mean.b_a = zeros(3*K,1);
init.cov.b_a = acc.init_b_a^2*eye(3*K)*3;

init.mean.b_g = zeros(3,1);
init.cov.b_g = eye(3)*gyro.init_b_g^2/K*3;

% Run Filters

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
simdata.accelerometer_array_2nd_order.r = pos.bodyCorr;
simdata.accelerometer_array_2nd_order.Q_acc = acc.noiseSigDisc^2*eye(3*K);
simdata.accelerometer_array_2nd_order.Q_bias_acc = acc.Q_b_a^2*eye(3*K);
simdata.accelerometer_array_2nd_order.R_pos = pos.sig^2*eye(3);
simdata.accelerometer_array_2nd_order.R_gyro = gyro.noiseSigDisc^2*eye(3)/K;
simdata.accelerometer_array_2nd_order.Q_bias_gyro = gyro.Q_b_g^2*eye(3)/K;
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
simdata.accelerometer_array_1st_order.r = pos.bodyCorr;
simdata.accelerometer_array_1st_order.Q_acc = acc.noiseSigDisc^2*eye(3*K);
simdata.accelerometer_array_1st_order.Q_bias_acc = acc.Q_b_a^2*eye(3*K);
simdata.accelerometer_array_1st_order.R_pos = pos.sig^2*eye(3);
simdata.accelerometer_array_1st_order.R_gyro = gyro.noiseSigDisc^2*eye(3);
simdata.accelerometer_array_1st_order.Q_bias_gyro = gyro.Q_b_g^2*eye(3)/K;
simdata.accelerometer_array_1st_order.do_gyro_updates = true;
simdata.accelerometer_array_1st_order.do_position_updates = true;
simdata.accelerometer_array_1st_order.do_rotation_updates = false;
simdata.accelerometer_array_1st_order.label = "1st order accelerometer array";
% 
sensorData = struct;
sensorData.acc_measurements = truth.acc_u + truth.biasAcc + acc.noise;
sensorData.gyro_measurements = truth.gyro_u + truth.biasGyro + gyro.noise;
sensorData.position_measurements = truth.pos + pos.noise;

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
[~,resSingle.accelerometer_array_1st_order] = run_filter(sensorData, init, simdata.accelerometer_array_1st_order, run_settings, S_true);
[~, resSingle.accelerometer_array_2nd_order] = run_filter(sensorData, init, simdata.accelerometer_array_2nd_order, run_settings, S_true);

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
%% plotResults
function plotResults(meas, resSingle)

plotBias = 0;
plotRotation = 0;
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


%%