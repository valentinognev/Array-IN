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
    % data=load('../../Data/lssimDatad10.mat');
    data=load('../../Data/genData.mat');
    [tim, world, body, npoints]=getBodyData(data);
    cardou12axis(body);
    
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
    gyro.noise = chol(gyro.Q)*randn(3,meas.N)*0;

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
    figure(234);
    plot(tim, squeeze(body.omega(1,1,:))); hold on; plot(tim, squeeze(body.omega(2,1,:))); plot(tim, squeeze(body.omega(3,1,:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END INJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[init, simdata, sensorData, run_settings] = setSimulationParameters(acc, gyro, pos, truth, settings_default);
[resSingle]=runSimulation(init, sensorData, simdata, run_settings, truth, S_true, w_b);
plotResults(meas, resSingle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getBodyData
function [tim, world, body, npoints]=getBodyData(dataS)
% cg.x=median(data.coordx(:,:), 2);
% cg.y=cg.x*0;
% cg.z=median(data.coordz(:,:), 2);
indz=dataS.tim==0;indz(1)=false;
dataS.disx(indz,:)=[];
dataS.disy(indz,:)=[];
dataS.disz(indz,:)=[];
dataS.velx(indz,:)=[];
dataS.vely(indz,:)=[];
dataS.velz(indz,:)=[];
dataS.accx(indz,:)=[];
dataS.accy(indz,:)=[];
dataS.accz(indz,:)=[];
dataS.coordx(indz,:)=[];
dataS.coordy(indz,:)=[];
dataS.coordz(indz,:)=[];
dataS.tim(indz)=[];
downsampling = 1;%1000;
indzd = 1:downsampling:length(dataS.tim);

timfactor = 1;%1e-6;
posfactor = 1;%1e-3;
velfactor = posfactor/timfactor;
accfactor = posfactor/timfactor/timfactor;

%  Downsampling
data.disxS=dataS.disx(indzd,:)'*posfactor; data.disyS=dataS.disy(indzd,:)'*posfactor; data.diszS=dataS.disz(indzd,:)'*posfactor;
data.velxS=dataS.velx(indzd,:)'*velfactor; data.velyS=dataS.vely(indzd,:)'*velfactor; data.velzS=dataS.velz(indzd,:)'*velfactor;
data.accxS=dataS.accx(indzd,:)'*accfactor; data.accyS=dataS.accy(indzd,:)'*accfactor; data.acczS=dataS.accz(indzd,:)'*accfactor;
data.coordx=dataS.coordx(indzd,:)'*posfactor; data.coordy=dataS.coordy(indzd,:)'*posfactor; data.coordz=dataS.coordz(indzd,:)'*posfactor;
data.tim=dataS.tim(indzd)*timfactor;
tim = data.tim;

dt=diff(tim'); dt(end+1)=dt(end);

% recalculation of velocities and accelerations
data.velxD = diff(data.coordx,1,2)./diff(data.tim');    data.velyD = diff(data.coordy,1,2)./diff(data.tim)';   data.velzD = diff(data.coordz,1,2)./diff(data.tim');
data.velxD = [data.velxD, data.velxD(:,end)];      data.velyD = [data.velyD, data.velyD(:,end)];       data.velzD = [data.velzD, data.velzD(:,end)];

data.accxD = diff(data.velxD,1,2)./diff(data.tim');    data.accyD = diff(data.velyD,1,2)./diff(data.tim');   data.acczD = diff(data.velzD,1,2)./diff(data.tim');
data.accxD = [data.accxD, data.accxD(:,end)];          data.accyD = [data.accyD, data.accyD(:,end)];         data.acczD = [data.acczD, data.acczD(:,end)];
data.accxSD = diff(data.velxS,1,2)./diff(data.tim');  data.accySD = diff(data.velyS,1,2)./diff(data.tim'); data.acczSD = diff(data.velzS,1,2)./diff(data.tim');
data.accxSD = [data.accxSD, data.accxSD(:,end)];       data.accySD = [data.accySD, data.accySD(:,end)];      data.acczSD = [data.acczSD, data.acczSD(:,end)];

data.disxD = cumsum(data.velxD.*dt,2);    data.disy = cumsum(data.velyD.*dt,2);    data.disz = cumsum(data.velzD.*dt,2);
data.velxDC = cumsum(data.accxD.*dt,2);   data.velyDC = cumsum(data.accyD.*dt,2);   data.velzDC = cumsum(data.acczD.*dt,2);
data.velxDC = data.velxDC+data.velxD(:,1);data.velyDC = data.velyDC+data.velyD(:,1);data.velzDC = data.velzDC+data.velzD(:,1);
data.velxSC = cumsum(data.accxS.*dt,2); data.velySC = cumsum(data.accyS.*dt,2); data.velzSC = cumsum(data.acczS.*dt,2);
data.velxSC = data.velxSC+data.velxD(:,1);data.velySC = data.velySC+data.velyD(:,1);data.velzSC = data.velzSC+data.velzD(:,1);

indp=5;
% plot(tim,data.velxD(indp,:)); hold on; plot(tim,data.velxS(indp,:)); plot(tim,data.velxDC(indp,:)); plot(tim,data.velxSC(indp,:)); 
% plot(tim,data.velzD(indp,:)); hold on; plot(tim,data.velzS(indp,:)); plot(tim,data.velzDC(indp,:)); plot(tim,-data.velzSC(indp,:)/10); 
% plot(tim,data.accxD(indp,:)); hold on; plot(tim,data.accxS(indp,:)); plot(tim,data.accxSD(indp,:)); plot(tim,data.accxDC(indp,:)); plot(tim,data.accxSC(indp,:)); 
% plot(tim,data.acczD(indp,:)); hold on; plot(tim,data.acczS(indp,:)); plot(tim,data.acczSD(indp,:)); plot(tim,data.acczDC(indp,:)); plot(tim,-data.acczSC(indp,:)/10); 
% plot(tim,data.accx(indp,:)); hold on; plot(tim,data.accxS(indp,:)); 
% plot(tim,data.disy(indp,:)); hold on; plot(tim,data.disyS(indp,:)); 
% plot(tim,data.disz(indp,:)); hold on; plot(tim,data.diszS(indp,:)); 
% legend({'D','S','DC','SC'})

% calculation body frame origin position and directions
% calculation of frame direction
npoints = size(data.coordx,1);
ntim = length(data.tim);

originNID = 1;%16;
xdirNID = 6;%15;
zdirNID = 2;%24;
dirX=[data.coordx(xdirNID,:)-data.coordx(originNID,:); 
      data.coordy(xdirNID,:)-data.coordy(originNID,:); 
      data.coordz(xdirNID,:)-data.coordz(originNID,:)];
d24_16=[data.coordx(zdirNID,:)-data.coordx(originNID,:); 
        data.coordy(zdirNID,:)-data.coordy(originNID,:); 
        data.coordz(zdirNID,:)-data.coordz(originNID,:)];
for i=1:ntim
    dirZ(:,i)=cross(dirX(:,i), cross(d24_16(:,i), dirX(:,i)));
    dirY(:,i)=cross(dirZ(:,i), dirX(:,i));
    dirX(:,i)=dirX(:,i)/norm(dirX(:,i));
    dirY(:,i)=dirY(:,i)/norm(dirY(:,i));
    dirZ(:,i)=dirZ(:,i)/norm(dirZ(:,i));
    disp ''
end
world.frame.x=dirX;  % direction of body frame X coordinate in world frame
world.frame.y=dirY;  % direction of body frame Y coordinate in world frame
world.frame.z=dirZ;  % direction of body frame Z coordinate in world frame

world.origin.pos = [data.coordx(originNID,:); data.coordy(originNID,:); data.coordz(originNID,:)];
world.origin.vel = [data.velxD(originNID,:); data.velyD(originNID,:); data.velzD(originNID,:)];
world.origin.acc = [data.accxSD(originNID,:); data.accySD(originNID,:); data.acczSD(originNID,:)];

% calculation of rotation matrices and quaternions : BW and WB 
world.rot_wb = zeros(3,3,ntim);
world.rot_bw = zeros(3,3,ntim);
world.quat_wb = zeros(4,ntim);
world.quat_bw = zeros(4,ntim);
for i=1:ntim
    world.rot_wb(:,:,i) = [dirX(:,i), dirY(:,i), dirZ(:,i)];
    world.rot_bw(:,:,i) = [dirX(:,i), dirY(:,i), dirZ(:,i)]';
    world.quat_wb(:,i) = dcm2quat(squeeze(world.rot_wb(:,:,i)));
    world.quat_bw(:,i) = dcm2quat(squeeze(world.rot_bw(:,:,i)));
    disp ''
end

% calculation of body angular velocity - omega
body.omega = angular_velocitiesMat(data.tim, world.rot_wb)';
body.omega_qwb = angular_velocities(data.tim, world.quat_wb);
body.omega_qwb2 = angular_velocities2(data.tim, world.quat_wb,'w');
body.omegadot_qwb = angular_accelerations(data.tim, world.quat_wb,'w');
body.omegadot = diff(body.omega,1,2); body.omegadot(:,end+1) = body.omegadot(:,end);
body.omegadot = body.omegadot./([1;1;1]*dt);

for i=1:ntim
    world.omega(:,i) = world.rot_wb(:,:,i)*body.omega(:,i);
    world.omegadot(:,i) = world.rot_wb(:,:,i)*body.omegadot(:,i); 
end

% Test of omegaB
% skewMat = @(omega) [0, -omega(3), omega(2);
%                      omega(3), 0, -omega(1);
%                     -omega(2), omega(1), 0];
% for i=1:(ntim-1)
%     R1=world.rot_wb(:,:,i);
%     R2=world.rot_wb(:,:,i+1);
%     omegaWmat=skewMat(world.omega(:,i));
%     deltaT=tim(i+1)-tim(i);
%     Rtest=expm(omegaWmat*deltaT)*R1-R2;
%     test(i)=norm(Rtest,"fro")/norm(R2,"fro")*100;
%     disp ''
% end


% plot(data.tim, world.omega_qwb(1,:));hold on;plot(data.tim, world.omega_qwb2(1,:));plot(data.tim, world.omegaMat(1,:));
% plot(data.tim, world.omega_qwb(2,:));hold on;plot(data.tim, world.omega_qwb2(2,:));plot(data.tim, world.omegaMat(2,:));
% plot(data.tim, world.omega_qwb(3,:));hold on;plot(data.tim, world.omega_qwb2(3,:));plot(data.tim, world.omegaMat(3,:));
%
% plot(data.tim, world.omegadot_qwb(2,:));plot(data.tim, world.omegadot_qwb(3,:));
% plot(data.tim, world.omegadot_2(2,:));plot(data.tim, world.omegadot_2(3,:));

body.pos = zeros(3,npoints, length(tim)); body.vel = zeros(3,npoints, length(tim)); body.acc = zeros(3,npoints, length(tim));
body.origin.pos = zeros(3, length(tim));  body.origin.vel = zeros(3, length(tim));  body.origin.acc = zeros(3, length(tim));
                                          body.origin.velS = zeros(3, length(tim)); body.origin.accS = zeros(3, length(tim));
world.pos = zeros(3,npoints, length(tim));
testacc=zeros(3, length(tim));
for i=1:ntim
    rot_bw = squeeze(world.rot_bw(:,:,i));
    world.pos(:,:,i)=[data.coordx(:,i)-world.origin.pos(1,i),...
                      data.coordy(:,i)-world.origin.pos(2,i),...
                      data.coordz(:,i)-world.origin.pos(3,i)]';
    world.vel(:,:,i)=[data.velxD(:,i), data.velyD(:,i), data.velzD(:,i)]';
    world.acc(:,:,i)=[data.accxSD(:,i), data.accySD(:,i), data.acczSD(:,i)]';

    body.pos(:,:,i) = rot_bw*world.pos(:,:,i);              
    body.vel(:,:,i) = rot_bw*world.vel(:,:,i);              
    body.acc(:,:,i) = rot_bw*world.acc(:,:,i);
    body.origin.pos(:,i) = rot_bw*world.origin.pos(:,i);    
    body.origin.vel(:,i) = rot_bw*world.origin.vel(:,i);    
    body.origin.acc(:,i) = rot_bw*world.origin.acc(:,i);

    % ind=10;
    % frame=world;
    % rb_ = frame.pos(:,ind,i);
    % omb = frame.omega(:,i);
    % omdotb = frame.omegadot(:,i);
    % sb = frame.origin.acc(:,i);
    % vb = frame.origin.vel(:,i);
    % fb_ = frame.acc(:,ind,i);
    % testvel(:,i) = vb + cross(omb, rb_);
    % testacc(:,i) = sb + cross(omdotb, rb_) + cross(omb, cross(omb, rb_));
    disp ''
end
% plot(tim,testvel(1,:)); hold on; plot(tim,squeeze(frame.vel(1,ind,:))); 
% plot(tim,testvel(3,:)); hold on; plot(tim,squeeze(frame.vel(3,ind,:))); 
% plot(tim,testvel(2,:)); hold on; plot(tim,squeeze(frame.vel(2,ind,:)));

% plot(tim,testacc(1,:)); hold on; plot(tim,squeeze(frame.acc(1,ind,:))); 
% plot(tim,testacc(3,:)); hold on; plot(tim,squeeze(frame.acc(3,ind,:))); 
% plot(tim,testacc(2,:)); hold on; plot(tim,squeeze(frame.acc(2,ind,:))); 

dpos = diff(world.pos,1,3);  dpos(:,:,end+1)=dpos(:,:,end);
ddpos = diff(dpos,1,3);  ddpos(:,:,end+1)=ddpos(:,:,end);
for i=1:ntim
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
[~, resSingle.accelerometer_array_1st_order] = run_filter(sensorData, init, simdata.accelerometer_array_1st_order, run_settings, S_true);
% [~, resSingle.accelerometer_array_2nd_order] = run_filter(sensorData, init, simdata.accelerometer_array_2nd_order, run_settings, S_true);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% cardou12axis
function cardou12axis(data)
fbarr = data.acc;                   % acceleration measurements
sbarr = data.origin.acc;                 % acceleration of body c.g.
ombarr = data.omega;                % omega body
ombdotarr = data.omegadot;          % omega dot body
rb = squeeze(data.pos(:,:,1));   % accelerometer position

sz=size(fbarr,3);
npoints=size(rb,2);
omb=ombarr(:,1);
for j=1:sz
    fb=squeeze(fbarr(:,:,j));
    omdotb=ombdotarr(:,j);
    sb=sbarr(:,j);

    for i = 1:npoints
        rb_ = rb(:,i);
        fb_ = fb(:,i);
        test(:,i) = sb + cross(omdotb, rb_) + cross(omb, cross(omb, rb_)) - fb_;
        percent(i)=norm(test(:,i))/norm(fb_)*100;
        % fb((i*3-2):(i*3)) = sb + cross(omdotb, rb_) + cross(omb, cross(omb, rb_));
    end

    acc = fb(:);
    evecbase = [ 1, 0, 0;
                 0, 1, 0;
                 0, 0, 1];
    evec = repmat(evecbase,npoints,1);

    mdofs = 3*npoints;
    rvec=zeros(3,mdofs);
    for i=1:npoints
        mat=[rb(1:3,i), rb(1:3,i), rb(1:3,i)];
        rvec(:,(3*i-2):(3*i))=mat;
    end

    CPM = @(x) [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
    Sigmaa = @(x) [0, -x(1), -x(1), 0, x(3), x(2); -x(2), 0, -x(2), x(3), 0, x(1); -x(3), -x(3), 0, x(2), x(1), 0];

    Ap = evec;
    R = zeros(3, 3*mdofs);
    F = zeros(mdofs, 3*mdofs);
    Sigma = zeros(6, 3*mdofs);

    for i = 1:mdofs
        try
        F(i,(3*i-2):(3*i)) = evec(i,:);
        R(:,(3*i-2):(3*i)) = CPM(rvec(:,i));
        Sigma(:,(3*i-2):(3*i)) = Sigmaa(rvec(:,i)).';
        catch
            disp ''
        end
    end

    At = F * R.';
    Ar = F * Sigma.';
    A = [Ap, At, Ar];

    res = lsqminnorm(A, acc);

    bdotdot = res(1:3);
    wdot = res(4:6);
    w0sq = res(7);
    w1sq = res(8);
    w2sq = res(9);
    w1w2 = res(10);
    w0w2 = res(11);
    w0w1 = res(12);

    Ws = [-w1sq-w2sq, w0w1, w0w2; w0w1, -w0sq-w2sq, w1w2; w0w2, w1w2, -w0sq-w1sq];

    wCANP = calcCANP(Ws, omb);
    wCAD = calcCAD(Ws, omb);
    wCAAD = calcCAAD(Ws, omb);
    wCAAM = calcCAAM(Ws, omb);
    omb=wCANP';
end

disp ''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% calcCANP
function wCANP = calcCANP(Ws, wTA)
if size(wTA,2)==3
    wTA=wTA';
end
mu2 = -(Ws(1,1)+Ws(2,2)+Ws(3,3));
s12 = Ws(1,2)^2;
s23 = Ws(2,3)^2;
s31 = Ws(3,1)^2;
s13 = Ws(1,2)*Ws(2,3);
d12 = Ws(1,1)*Ws(2,2);
mu1 = d12 + Ws(2,2)*Ws(3,3) + Ws(1,1)*Ws(3,3) - s12 - s23 - s31;
mu0 = -d12*Ws(3,3) - 2*s13*Ws(3,1) + s12*Ws(3,3) + s23*Ws(1,1) + s31*Ws(2,2);
ni2 = mu2/3;
theta2 = ni2^2;
q = mu1/3 - theta2;

if q >= 0
    wCANP = [0, 0, 0];
    return;
else
    r = (mu1*ni2 - mu0) / 2 - ni2*theta2;
    alpha = sqrt(-q);
    beta = alpha^3;
end

if beta <= r
    wCANP = [0, 0, 0];
    return;
else
    lam = 2 * alpha * cos(acos(r / beta) / 3) - ni2;
    delta = (lam + mu2) / 2;
end

if delta <= 0 || (lam * mu0) > 0
    wCANP = [0, 0, 0];
    return;
else
    wcanpnorm = sqrt(delta);
    zeta11 = Ws(1,1) - lam;
    zeta22 = Ws(2,2) - lam;
    zeta33 = Ws(3,3) - lam;
    ksi11 = zeta22 * zeta33 - s23;
    ksi22 = zeta33 * zeta11 - s31;
    ksi33 = zeta11 * zeta22 - s12;
    ksi12 = Ws(2,3) * Ws(3,1) - Ws(1,2) * Ws(3,3);
    ksi23 = Ws(1,2) * Ws(3,1) - Ws(2,3) * Ws(1,1);
    ksi31 = s13 - Ws(3,1) * Ws(2,2);
    adjX = [ksi11, ksi12, ksi31; ksi12, ksi22, ksi23; ksi31, ksi23, ksi33];
    v = adjX * wTA;
    vunit = v / norm(v);
    wCANP = wcanpnorm * vunit.';
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% calcCAD
function wCAD = calcCAD(Ws, wTA)
trWs = trace(Ws);
if trWs < 0 && sum(abs(wTA)) > 0
    zeta0 = Ws(1,1) - 0.5 * trWs;
    zeta1 = Ws(2,2) - 0.5 * trWs;
    zeta2 = Ws(3,3) - 0.5 * trWs;
    wCAD0 = sign(wTA(1)) * heaviside(zeta0) * sqrt(zeta0);
    wCAD1 = sign(wTA(2)) * heaviside(zeta1) * sqrt(zeta1);
    wCAD2 = sign(wTA(3)) * heaviside(zeta2) * sqrt(zeta2);
    wCAD = [wCAD0, wCAD1, wCAD2];
    return;
else
    wCAD = [0, 0, 0];
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% calcCAAD
function wCAAD = calcCAAD(Ws, wTA)
if size(wTA,2)==3
    wTA=wTA';
end
trWs = trace(Ws);
if trWs < 0 && sum(abs(wTA)) > 0
    adjWs = inv(Ws) * det(Ws);
    wCAADnorm = sqrt(-0.5 * trWs);
    v = adjWs * wTA;
    vunit = v / norm(v);
    wCAAD = wCAADnorm * vunit.';
    return;
else
    wCAAD = [0, 0, 0];
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% calcCAAM
function wCAAM = calcCAAM(Ws, wTA)
if size(wTA,1)==3
    wTA=wTA';
end
trWs = trace(Ws);
if trWs < 0
    adjWs = inv(Ws) * det(Ws);
    wCAAMnorm = sqrt(-0.5 * trWs);
    Xtop = Ws;
    Xbot = wTA * adjWs;
    X = [Xtop; Xbot];
    Y = [0, 0, 0, (-trWs/2)^3];
    [Q, R] = qr(X);
    uw = R \ (Q.' * Y.');
    uwunit = uw / norm(uw);
    wCAAM = wCAAMnorm * uwunit.';
    return;
else
    wCAAM = [0, 0, 0];
    return;
end


