function test
% R1 = [ 0 1 0; 1 0 0; 0 0 -1];
% R2 = [ 0.4330    0.2500   -0.8660; ...
%        0.1768    0.9186    0.3536; ...
%        0.8839   -0.3062    0.3536];
R1=[ 1     0     0;
     0     1     0;
     0     0     1];

R2=[                     1                         0                         0;
                         0            0.999992000016      -0.00399998800002001;
                         0       0.00399998800002001            0.999992000016];

q(1) = quaternion(dcm2quat(R1));
q(2) = quaternion(dcm2quat(R2));

v=[1;2;3];  v=v/norm(v);

v2_dcm_1 = R1*v;
v2_dcm_2 = R2*v;

qv = quaternion([0;v]');
v2_q_1 = quatconj(q(1))*qv*q(1);
v2_q_2 = quatconj(q(2))*qv*q(2);
qx=quatconj(q(1))*q(2);
v2_q_2_2=quatconj(qx)*v2_q_1*qx;

deltamat=quat2dcm(qx);
testdelta = deltamat*R1*v;

dR = R2-R1;
dt=.001;
i=1;
omegaMat = (dR/dt)*R1';
omegaMat = (omegaMat-omegaMat')/2;
omegaRot = [omegaMat(3,2), omegaMat(1,3), omegaMat(2,1)];
[omegaQ] = angular_velocities2(dt, q(1), q(2), 'w');
skew = @(omega) [0, -omega(3), omega(2);
                 omega(3), 0, -omega(1);
                -omega(2), omega(1), 0];

RtestMat=expm(skew(omegaRot)*dt)*R1;
RtestQ=expm(skew(omegaRot)*dt)*R1;

disp ''

function [omega] = angular_velocities2(dt, q1, q2, frame)
% frame: 'w' for world frame, 'b' for body frame
if nargin<3
    frame = 'w';
end

qdot = (q2-q1)*(1/dt);  
[q_1,q_2,q_3,q_4]=parts(qdot);
qdot=[q_1,q_2,q_3,q_4]';

if frame=='w'
    rotmat = QmatbarT(q1);
else
    rotmat = QmatT(q1);
end
omega = 2*rotmat*qdot;  %2*QmatT(qHist)*qdot';
omega = sign(omega(1))*omega(2:4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Qmat
function [qmat] = Qmat(q)
qmat = zeros(4,4);
[q1,q2,q3,q4]=parts(q);
q=[q1,q2,q3,q4];
qmat(1,1) = q(1);  qmat(1,2) = -q(2); qmat(1,3) = -q(3); qmat(1,4) = -q(4);
qmat(2,1) = q(2);  qmat(2,2) =  q(1); qmat(2,3) = -q(4); qmat(2,4) =  q(3);
qmat(3,1) = q(3);  qmat(3,2) =  q(4); qmat(3,3) =  q(1); qmat(3,4) = -q(2);
qmat(4,1) = q(4);  qmat(4,2) = -q(3); qmat(4,3) =  q(2); qmat(4,4) =  q(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QmatT
function [qmatT] = QmatT(q)
qmatT = zeros(4,4);
[q1,q2,q3,q4]=parts(q);
q=[q1,q2,q3,q4];
qmatT(1,1) =  q(1);  qmatT(1,2) =  q(2); qmatT(1,3) =  q(3); qmatT(1,4) =  q(4);
qmatT(2,1) = -q(2);  qmatT(2,2) =  q(1); qmatT(2,3) =  q(4); qmatT(2,4) = -q(3);
qmatT(3,1) = -q(3);  qmatT(3,2) = -q(4); qmatT(3,3) =  q(1); qmatT(3,4) =  q(2);
qmatT(4,1) = -q(4);  qmatT(4,2) =  q(3); qmatT(4,3) = -q(2); qmatT(4,4) =  q(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Qmatbar
function [qmatbar] = Qmatbar(q)
qmatbar = zeros(4,4);
[q1,q2,q3,q4]=parts(q);
q=[q1,q2,q3,q4];
qmatbar(1,1) = q(1);  qmatbar(1,2) = -q(2); qmatbar(1,3) = -q(3); qmatbar(1,4) = -q(4);
qmatbar(2,1) = q(2);  qmatbar(2,2) =  q(1); qmatbar(2,3) =  q(4); qmatbar(2,4) = -q(3);
qmatbar(3,1) = q(3);  qmatbar(3,2) = -q(4); qmatbar(3,3) =  q(1); qmatbar(3,4) =  q(2);
qmatbar(4,1) = q(4);  qmatbar(4,2) =  q(3); qmatbar(4,3) = -q(2); qmatbar(4,4) =  q(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QmatbarT
function [qmatbarT] = QmatbarT(q)
qmatbarT = zeros(4,4);
[q1,q2,q3,q4]=parts(q);
q=[q1,q2,q3,q4];
qmatbarT(1,1) =  q(1);  qmatbarT(1,2) =  q(2); qmatbarT(1,3) =  q(3); qmatbarT(1,4) =  q(4);
qmatbarT(2,1) = -q(2);  qmatbarT(2,2) =  q(1); qmatbarT(2,3) = -q(4); qmatbarT(2,4) =  q(3);
qmatbarT(3,1) = -q(3);  qmatbarT(3,2) =  q(4); qmatbarT(3,3) =  q(1); qmatbarT(3,4) = -q(2);
qmatbarT(4,1) = -q(4);  qmatbarT(4,2) = -q(3); qmatbarT(4,3) =  q(2); qmatbarT(4,4) =  q(1);
