function test
% R1 = [ 0 1 0; 1 0 0; 0 0 -1];
% R2 = [ 0.4330    0.2500   -0.8660; ...
%        0.1768    0.9186    0.3536; ...
%        0.8839   -0.3062    0.3536];
R1=[         0.999999935000001     -0.000299989993000128      0.000200014995333159
      0.000300009992999895         0.999999950000001     -9.99699976670376e-05
     -0.000199984995333509      0.000100029997666338               0.999999975];

R2=[   0.999999740000012     -0.000599959944001887      0.000400059962663867
      0.000600039943998154         0.999999800000009     -0.000199879981338913
     -0.000399939962669467      0.000200119981327713         0.999999900000005];

q(1) = quaternion(dcm2quat(R1));
q(2) = quaternion(dcm2quat(R2));

v=[1;2;3];  v=v/norm(v);

v2_dcm_1 = R1*v;
v2_dcm_2 = R2*v;

qv = quaternion([0;v]');
v2_q_1 = quatconj(q(1))*qv*q(1);
v2_q_2 = quatconj(q(2))*qv*q(2);

dR = R2-R1;
dt=1;
i=1;
omegaMat = (dR/dt)*R1';
omegaMat = (omegaMat-omegaMat')/2;
omegaRot = [omegaMat(3,2), omegaMat(1,3), omegaMat(2,1)];

Rtest=expm(omegaMat*dt)*R1;

disp ''

function [omega] = angular_velocities2(dt, q1, q2, frame)
% frame: 'w' for world frame, 'b' for body frame
if nargin<3
    frame = 'w';
end

qdot = (q2-q1)/dt;  

if frame=='w'
    rotmat = QmatbarT(qHist);
else
    rotmat = QmatT(qHist);
end

omega = 2*rotmat*qdot;  %2*QmatT(qHist)*qdot';
omega = omega(2:4,:);
