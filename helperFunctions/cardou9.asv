function [bdotdot, omega_dot_b, omega_b]=cardou9(acc_b, omega_b_est, pos_b, dt)
% Cardou, 2010, Computing the Rigid body acceleration field from nine accelerometer measurements

% IN:
%    acc_b - acceleration measurements
%    ombest - omega body estimation
% OUT:
%    bdotdot - acceleration of body c.g.
%    omega_dot_b - angular acceleration
%    omega_b - angular velocity
bdotdot = [0.9981, -66.6620, 1.0000];
omega_dot_b = [2.9999, 2.0002, 0.9999];
omega = [0.1022, 0.2015, 0.3007];

pos_b=pos_b(:,[1,2,4]);
acc_b=acc_b([1:6,10:12]);
% bdotdot=[]; omega_dot_b=[]; omega_b=[];
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
    F(i,(3*i-2):(3*i)) = evec(i,:);
    R(:,(3*i-2):(3*i)) = CPM(rvec(:,i));
    Sigma(:,(3*i-2):(3*i)) = Sigmaa(rvec(:,i)).';
end

At = F * R.';
Ar = F * Sigma.';

% Commented out parts that were specific to Python and need adaptation
% % Ar = Ar.' * rho;
A = [Ap, At, Ar];
[Q, S] = qr(A);

Q1 = Q(:, 1:6);
Q2 = Q(:, 7:9);
if size(Q, 2) > 9
    Q3 = Q(:, 10:end);
else
    Q3 = [];
end

S11 = S(1:6, 1:6);
S22 = S(7:9, 7:12);
S12 = S(1:6, 7:12);
O3x6 = S(7:9, 1:6);
if size(S, 2) > 12
    S32 = S(:, 7:12);
    Onm9x6 = S(:, 1:6);
else
    S32 = [];
    Onm9x6 = [];
end

rankA = rank(A);
Q2a = Q2.' * acc;

s11 = S22(1, 1); s12 = S22(1, 2); s13 = S22(1, 3); s14 = S22(1, 4); s15 = S22(1, 5); s16 = S22(1, 6);
s22 = S22(2, 2); s23 = S22(2, 3); s24 = S22(2, 4); s25 = S22(2, 5); s26 = S22(2, 6);
s33 = S22(3, 3); s34 = S22(3, 4); s35 = S22(3, 5); s36 = S22(3, 6);
q7a = Q2a(1); q8a = Q2a(2); q9a = Q2a(3);

[v0, v1, v2, v3, v4] = getVcoeffs(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a);
zeta = roots([v4, v3, v2, v1, v0]);

[zeta_0, zeta_1, zeta_2, zeta_3] = quarticEquationRoots(v0, v1, v2, v3, v4);
zeta2=[zeta_0, zeta_1, zeta_2, zeta_3]';
w3arr = [];

for i = 1:length(zeta)
    w3p = sqrt(zeta(i));
    w3n = -sqrt(zeta(i));
    w3arr = [w3arr, w3p, w3n];
end

wArr = [];
wdotArr = [];
bddArr = [];

for w3 = w3arr
    [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10] = getuvector(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a, w3);

    C1 = [s11, s16, s15*w3, s12, s14*w3;
        0, s26, s25*w3, s22, s24*w3;
        0, s36, s35*w3, 0, s34*w3;
        3*u1, 2*u2, 2*u3, u4, u5;
        u2, 2*u4, u5, 3*u7, 2*u8;
        u3, u5, 2*u6, u8, 2*u9];

    c2 = -[s13*w3^2 - q7a;
        s23*w3^2 - q8a;
        s33*w3^2 - q9a;
        u6;
        u9;
        3*u10];

    rankC1 = rank(C1);
    fprintf('rankC1: %d\n', rankC1);

    [Q_2, S_2] = qr(C1);
    lhs = Q_2.' * c2;
    ww4 = lhs(5) / S_2(5, 5);
    ww3 = (lhs(4) - S_2(4, 5) * ww4) / S_2(4, 4);
    ww2 = (lhs(3) - S_2(3, 4) * ww3 - S_2(3, 5) * ww4) / S_2(3, 3);
    ww1 = (lhs(2) - S_2(2, 3) * ww2 - S_2(2, 4) * ww3 - S_2(2, 5) * ww4) / S_2(2, 2);
    ww0 = (lhs(1) - S_2(1, 2) * ww1 - S_2(1, 3) * ww2 - S_2(1, 4) * ww3 - S_2(1, 5) * ww4) / S_2(1, 1);
    wwi = [ww0, ww1, ww2, ww3, ww4];

    if ww0 < 0 || ww3 < 0
        continue;
    end

    w1abs = real(sqrt(ww0));
    w2abs = real(sqrt(ww3));

    if ismembertol(abs(real(ww2)), w1abs) && ismembertol(abs(real(ww4)), w2abs)
        wArr = [wArr; real([ww2, ww4, w3])];
    elseif real(ww1) >= 0
        wArr = [wArr; [w1abs, w2abs, real(w3)]; [-w1abs, -w2abs, real(w3)]];
    else
        wArr = [wArr; [w1abs, -w2abs, real(w3)]; [-w1abs, w2abs, real(w3)]];
    end
end

wArr = unique(real(wArr), 'rows');

for i = 1:size(wArr, 1)
    w1 = wArr(i, 1);
    w2 = wArr(i, 2);
    w3 = wArr(i, 3);

    ksii = [w1^2, w2^2, w3^2, w2*w3, w1*w3, w1*w2];
    Q1a = Q1.' * acc;
    S12ksii = S12 * ksii.';
    Q1amS12ksii = Q1a - S12ksii;

    xpt5 = Q1amS12ksii(6) / S11(6, 6);
    xpt4 = (Q1amS12ksii(5) - S11(5, 6) * xpt5) / S11(5, 5);
    xpt3 = (Q1amS12ksii(4) - S11(4, 5) * xpt4 - S11(4, 6) * xpt5) / S11(4, 4);
    xpt2 = (Q1amS12ksii(3) - S11(3, 4) * xpt3 - S11(3, 5) * xpt4 - S11(3, 6) * xpt5) / S11(3, 3);
    xpt1 = (Q1amS12ksii(2) - S11(2, 3) * xpt2 - S11(2, 4) * xpt3 - S11(2, 5) * xpt4 - S11(2, 6) * xpt5) / S11(2, 2);
    xpt0 = (Q1amS12ksii(1) - S11(1, 2) * xpt1 - S11(1, 3) * xpt2 - S11(1, 4) * xpt3 - S11(1, 5) * xpt4 - S11(1, 6) * xpt5) / S11(1, 1);
    xpti = [xpt0, xpt1, xpt2, xpt3, xpt4, xpt5];

    bdotdot = xpti(1:3);
    wdot = xpti(4:6);
    wdotArr = [wdotArr; wdot];
    bddArr = [bddArr; bdotdot];
end

wdotArr = real(wdotArr);
bddArr = real(bddArr);

for i = 1:size(wArr, 1)
    omdotb_ = wdotArr(i, :);
    omb_ = wArr(i, :);
    sb_ = bddArr(i, :);

    ksii = [omb_(1)^2, omb_(2)^2, omb_(3)^2, omb_(2) * omb_(3), omb_(1) * omb_(3), omb_(1) * omb_(2)];
    test = Ap * sb_.' + At * omdotb_.' + Ar * ksii.' - acc;
    % additional tests or operations can go here...
    test2 = [];

    for j = 1:length(rvec)
        rb_ = rvec(j, :);
        fb_ = acc(j);
        test2 = [test2; dot(evec(j, :), sb_ + cross(omdotb_, rb_) + cross(omb_, cross(omb_, rb_))) - fb_];
    end
end
disp ''

function roots = cubicEquationRoots(v0, v1, v2, v3)
zeta0 = -(v2/(3*v3)) - (2^(1/3)*(-v2^2 + 3*v1*v3))/(3*v3*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3)) + (1/(3*2^(1/3)*v3))*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3);
zeta1 = -(v2/(3*v3)) + ((1 + 1i*sqrt(3))*(-v2^2 + 3*v1*v3))/(3*2^(2/3)*v3*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3)) - (1/(6*2^(1/3)*v3))*((1 - 1i*sqrt(3))*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3));
zeta2 = -(v2/(3*v3)) + ((1 - 1i*sqrt(3))*(-v2^2 + 3*v1*v3))/(3*2^(2/3)*v3*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3)) - (1/(6*2^(1/3)*v3))*((1 + 1i*sqrt(3))*(-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2 + sqrt(4*(-v2^2 + 3*v1*v3)^3 + (-2*v2^3 + 9*v1*v2*v3 - 27*v0*v3^2)^2))^(1/3));
roots = [zeta0, zeta1, zeta2];

function [zeta0, zeta1, zeta2, zeta3] = quarticEquationRoots(v0, v1, v2, v3, v4)
sqrtV = sqrt(-4 * (v2^2 - 3 * v1 * v3 + 12 * v0 * v4)^3 + (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4)^2 + 0i);

sqrtA2_ = 0i + v3^2 / (2 * v4^2) - (4 * v2) / (3 * v4) - (2^(1/3) * (v2^2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3)) - (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3) / (3 * 2^(1/3) * v4);

sqrtA3_ = (-(v3^3 / v4^3) + (4 * v2 * v3) / v4^2 - (8 * v1) / v4) / (4 * sqrt(v3^2 / (4 * v4^2) - (2 * v2) / (3 * v4) + (2^(1/3) * (v2^2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3)) + (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3) / (3 * 2^(1/3) * v4)));

sqrtA1 = sqrt(0i + v3^2 / (4 * v4^2) - (2 * v2) / (3 * v4) + (2^(1/3) * (v2^2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3)) + (2 * v2^3 - 9 * v1 * v2 * v3 + 27 * v0 * v3^2 + 27 * v1^2 * v4 - 72 * v0 * v2 * v4 + sqrtV)^(1/3) / (3 * 2^(1/3) * v4));

zeta0 = -(v3 / (4 * v4)) - (1/2) * sqrtA1 - (1/2) * sqrt(sqrtA2_ - sqrtA3_);
zeta1 = -(v3 / (4 * v4)) - (1/2) * sqrtA1 + (1/2) * sqrt(sqrtA2_ - sqrtA3_);
zeta2 = -(v3 / (4 * v4)) + (1/2) * sqrtA1 - (1/2) * sqrt(sqrtA2_ + sqrtA3_);
zeta3 = -(v3 / (4 * v4)) + (1/2) * sqrtA1 + (1/2) * sqrt(sqrtA2_ + sqrtA3_);

function [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10] = getuvector(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a, w3)
u1 = (2*s11*s26*s35 - 2*s11*s25*s36)*w3;
u2 = (2*s11*s26*s34 + 4*s11*s22*s35 - 2*s11*s24*s36)*w3;
u3 = (4*s11*s26*s33 - 2*s11*s25*s34 + 2*s11*s24*s35 - 4*s11*s23*s36)*w3^2 + 4*q8a*s11*s36 - 4*q9a*s11*s26;
u4 = (4*s11*s22*s34 + 2*s16*s22*s35 - 2*s12*s26*s35 - 2*s15*s22*s36 + 2*s12*s25*s36)*w3;
u5 = (8*s11*s22*s33 - 2*s16*s25*s34 + 2*s15*s26*s34 + 2*s16*s24*s35 - 2*s14*s26*s35 - 2*s15*s24*s36 + 2*s14*s25*s36)*w3^2 - 8*q9a*s11*s22;
u6 = (4*s11*s24*s33 - 2*s16*s25*s33 + 2*s15*s26*s33 - 4*s11*s23*s34 + 2*s16*s23*s35 - 2*s13*s26*s35 - 2*s15*s23*s36 + 2*s13*s25*s36)*w3^3 + ...
    (-4*q9a*s11*s24 + 2*q9a*s16*s25 - 2*q9a*s15*s26 + 4*q8a*s11*s34 - 2*q8a*s16*s35 + 2*q7a*s26*s35 + 2*q8a*s15*s36 - 2*q7a*s25*s36)*w3;
u7 = (2*s16*s22*s34 - 2*s12*s26*s34 - 2*s14*s22*s36 + 2*s12*s24*s36)*w3;
u8 = (4*s16*s22*s33 - 4*s12*s26*s33 + 2*s15*s22*s34 - 2*s12*s25*s34 - 2*s14*s22*s35 + 2*s12*s24*s35 - 4*s13*s22*s36 + 4*s12*s23*s36)*w3^2 + ...
    (-4*q9a*s16*s22 + 4*q9a*s12*s26 - 4*q8a*s12*s36 + 4*q7a*s22*s36);
u9 = (4*s15*s22*s33 + 2*s16*s24*s33 - 4*s12*s25*s33 - 2*s14*s26*s33 - 2*s16*s23*s34 + 2*s13*s26*s34 - 4*s13*s22*s35 + 4*s12*s23*s35 + 2*s14*s23*s36 - 2*s13*s24*s36)*w3^3 + ...
    (-4*q9a*s15*s22 - 2*q9a*s16*s24 + 4*q9a*s12*s25 + 2*q9a*s14*s26 + 2*q8a*s16*s34 - 2*q7a*s26*s34 - 4*q8a*s12*s35 + 4*q7a*s22*s35 - 2*q8a*s14*s36 + 2*q7a*s24*s36)*w3;
u10 = -2*q9a*s15*s24*w3^2 + 2*q9a*s14*s25*w3^2 + 2*q8a*s15*s34*w3^2 - 2*q7a*s25*s34*w3^2 - 2*q8a*s14*s35*w3^2 + 2*q7a*s24*s35*w3^2 + 2*s15*s24*s33*w3^4 - ...
    2*s14*s25*s33*w3^4 - 2*s15*s23*s34*w3^4 + 2*s13*s25*s34*w3^4 + 2*s14*s23*s35*w3^4 - 2*s13*s24*s35*w3^4;


% acc=[16./175, 27./175, 4./175, 108./455, 48./455, 3./175, -6./35, -3./70, -19./70]'*pi^2;
% rho=0.1;
% evec=[0,1,0;
%       0,0,1;
%       1,0,0;
%       0,0,1;
%       0,1,0;
%       1,0,0;
%       1/sqrt(2),0,1/sqrt(2);
%       0,1/sqrt(2),1/sqrt(2);
%       1/sqrt(2),1/sqrt(2),0];
% 
% rvec=[rho,0,0;
%     rho,0,0;
%     0,rho,0;
%     0,rho,0;
%     0,0,rho;
%     0,0,rho;
%     rho/2,0,rho/2;
%     0,rho/2,rho/2;
%     rho/2,rho/2,0];
% 
% omb=[-4.4055, 1.6793, 1.5247, 3.7418, 4.4055, -1.6793, -1.5247, -3.7418;
%     -4.5368, 3.3585, 3.4270, 0.1900, 4.5368, -3.3585, -3.4270, -0.1900;
%     2.5476, 5.0378, 5.0563, 5.2394, -2.5476, -5.0378, -5.0563, -5.2394]';
% 
% omdotb=[-27.5229, 6.5074, 6.1716, 22.5814;
%     36.0568, -6.7677, -7.5911, 4.2272;
%     -5.4103, 3.3839, 4.5431, -23.6850]';
% 
% sb=[-2.3141, 0.0000, 0.1574,-2.2140;
%     -0.5553, -0.0000, -0.0745, 3.1998;
%     6.2508, 0.0000, -0.0073, -0.0150]';
