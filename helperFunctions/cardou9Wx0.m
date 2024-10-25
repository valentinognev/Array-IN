function [bdotdot, omega_dot_b, omega_b]=cardou9Wx0(acc_b, omega_b_est, pos_b, dt)
% IN:
%    acc_b - acceleration measurements
%    ombest - omega body estimation
% OUT:
%    bdotdot - acceleration of body c.g.
%    omega_dot_b - angular acceleration
%    omega_b - angular velocity
if size(pos_b,2)>3
    inds = [1,2,4];
    pos_b=pos_b(:,inds);
    acc_b = acc_b([1,2,3,4,5,6,10,11,12]);
end
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
A = [Ap, At, Ar];
A(:,[4,7,11,12]) = [];
res = lsqminnorm(A, acc);

bdotdot = res(1:3);
omega_dot_b = [0; res(4:5)];
w1sq = res(6);
w2sq = res(7);
w1w2 = res(8);

if w1sq<0 || w2sq<0
    omega_b = [0;0;0];
    return
end
w1=0;
if w1sq>0
    w1 = sqrt(w1sq);
end
w2=0;
if w2sq>0
    w2 = sqrt(w2sq);
end
if norm(omega_b_est)<1e8
    omega_b_est=omega_dot_b;
end
if w1w2>0
    omega_b = [0, w1, w2];
else
    omega_b = [0, w1, -w2];
end
if dot(omega_b, omega_b_est)<0
    omega_b = -omega_b;
end
omega_b = [0, sqrt(w1sq), sqrt(w2sq)]';
%%% test

ksi=[bdotdot', omega_dot_b(2),omega_dot_b(2), omega_b(2)^2, omega_b(3)^2, omega_b(2)*omega_b(3)];
test=(A*ksi'-acc)./acc*100;

disp ''
