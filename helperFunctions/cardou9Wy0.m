function [bdotdot, omega_dot_b, omega_b]=cardou9Wy0(acc_b, omega_b_est, pos_b, dt)
% IN:
%    acc_b - acceleration measurements
%    ombest - omega body estimation
% OUT:
%    bdotdot - acceleration of body c.g.
%    omega_dot_b - angular acceleration
%    omega_b - angular velocity

rb = pos_b';   % accelerometer position
acc = acc_b(1:9);

A =  [1, 0, 0,    0,    -rb(1,2),   0,     rb(1,3), -rb(1,1);
      0, 1, 0, -rb(1,3), rb(1,1),-rb(1,2),   0    , -rb(1,2);
      0, 0, 1,  rb(1,2),   0,    -rb(1,3), rb(1,1),     0   ;
      1, 0, 0,    0,    -rb(2,2),   0,     rb(2,3), -rb(2,1);
      0, 1, 0, -rb(2,3), rb(2,1),-rb(2,2),   0    , -rb(2,2);
      0, 0, 1,  rb(2,2),   0,    -rb(2,3), rb(2,1),     0   ;
      1, 0, 0,    0,    -rb(3,2),   0,     rb(3,3), -rb(3,1);
      0, 1, 0, -rb(3,3), rb(3,1),-rb(3,2),   0    , -rb(3,2);
      0, 0, 1,  rb(3,2),   0,    -rb(3,3), rb(3,1),     0   ];

res = lsqminnorm(A, acc);

bdotdot = res(1:3);
omega_dot_b = [res(4), 0, res(5)];
w0sq = res(6);
w0w2 = res(7);
w2sq = res(8);

w0=0;
if w0sq>0
    w0 = sqrt(w0sq);
end
w2=0;
if w2sq>0
    w2 = sqrt(w2sq);
end
if norm(omega_b_est)<1e8
    omega_b_est=omega_dot_b;
end
if w0w2>0
    omega_b = [w0, 0, w2];
else
    omega_b = [w0, 0, -w2];
end
if dot(omega_b(1,:),omega_b_est)<0
    omega_b = -omega_b;
end

disp ''
