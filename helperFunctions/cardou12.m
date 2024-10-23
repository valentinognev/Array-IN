function [bdotdot, omega_dot_b, omega_b]=cardou12(acc_b, omega_b_est, pos_b, dt)
% IN:
%    acc_b - acceleration measurements
%    ombest - omega body estimation
% OUT:
%    bdotdot - acceleration of body c.g.
%    omega_dot_b - angular acceleration
%    omega_b - angular velocity

rb = pos_b;   % accelerometer position

npoints=size(rb,2);

acc = acc_b(:);
evecbase = [1, 0, 0; 
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
omega_dot_b = res(4:6);
w0sq = res(7);
w1sq = res(8);
w2sq = res(9);
w1w2 = res(10);
w0w2 = res(11);
w0w1 = res(12);

Ws = [-w1sq-w2sq, w0w1, w0w2; w0w1, -w0sq-w2sq, w1w2; w0w2, w1w2, -w0sq-w1sq];

omb = omega_b_est;
if norm(omega_b_est)<1e-6
    omb=omega_dot_b*dt;
end

wCANP = calcCANP(Ws, omb)';
wCAD = calcCAD(Ws, omb)';
wCAAD = calcCAAD(Ws, omb)';
wCAAM = calcCAAM(Ws, omb)';
omega_b = wCANP;


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
