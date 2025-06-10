function [f,J] = dAdr(a,r,w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N_a = size(r,2); % Number of acc triads


R_skew = cell(N_a,1);
for k = 1:N_a
    R_skew{k} = skew_sym(r(:,k));
end
H = [-cat(1,R_skew{:}), repmat(eye(3),N_a,1)];


O2 = skew_sym(w)^2;
h = reshape(O2*r, [],1);
a_1 = a - h;
A_1 = H'*H;
A_1_inv = inv(A_1);

f = A_1\H'*a_1;

assert(length(f) == 6)
A_2 = kron(eye(N_a),O2);
A_3 = kron(f', -A_1_inv);
A_4 = kron(a_1', A_1_inv);
A_5 = -A_1\H'*A_2;

K = commutation_matrix(3*N_a,6);
A_6 = kron(H', eye(6))*K;
A_7 = kron(eye(6), H');
A_8 = A_3*A_6 + A_3*A_7 + A_4*K;

B_1 = skew_sym([-1 0 0]);
B_2 = skew_sym([0 -1 0]);
B_3 = skew_sym([0 0 -1]);

A_9 = [kron(eye(N_a), -B_1);
    kron(eye(N_a), -B_2);
    kron(eye(N_a), -B_3);
    zeros(9*N_a, 3*N_a)];

J = A_8 * A_9 + A_5;

assert(all(size(J) == [6,3*N_a]))

end

























function [d_h_d_omega] = d_h_d_omega(omega,r)
%D_H_D_OMEGA Summary of this function goes here
%   Detailed explanation goes here
N_a = size(r,2);
d_h_d_omega_parts = cell(N_a,1);
omega_hat = HatSO3(omega);
for k = 1:N_a
    r_k = r(:,k);
    d_h_d_omega_parts{k} = (-HatSO3(omega_hat*r_k) - omega_hat*HatSO3(r_k));
end
d_h_d_omega = cat(1, d_h_d_omega_parts{:});

end

function [d_h_d_omega] = d_h_d_omega_opt(w,r)
%D_H_D_OMEGA Summary of this function goes here
%   Detailed explanation goes here
N_a = size(r,2);
row1 = 1:3:3*N_a;
row2 = 2:3:3*N_a;
row3 = 3:3:3*N_a;
d_h_d_omega = zeros(3*N_a,3);
r1 = r(1,:);
r2 = r(2,:);
r3 = r(3,:);

r1w1 = w(1).*r1;
r1w2 = r1.*w(2);
r1w3 = r1.*w(3);

r2w1 = r2.*w(1);
r2w2 = w(2).*r2;
r2w3 = r2.*w(3);

r3w1 = r3.*w(1);
r3w2 = r3.*w(2);
r3w3 = w(3).*r3;

d_h_d_omega(row1,1) = r2w2 + r3w3;
d_h_d_omega(row2,1) = r1w2 - 2*r2w1;
d_h_d_omega(row3,1) = r1w3 - 2*r3w1;

d_h_d_omega(row1,2) = r2w1 - 2*r1w2;
d_h_d_omega(row2,2) = r1w1 + r3w3;
d_h_d_omega(row3,2) = r2w3 - 2*r3w2;

d_h_d_omega(row1,3) = r3w1 - 2*r1w3;
d_h_d_omega(row2,3) = r3w2 - 2*r2w3;
d_h_d_omega(row3,3) = r1w1 + r2w2;


end

function g = get_norm_g_kth()

g = 9.8183037; %  Lantmateriet, m/s^2

endfunction [y] = get_triad_form(x)
%GET_TRIAD_FORM Summary of this function goes here
%   Detailed explanation goes here
N_sens = size(x,1);
assert(mod(N_sens,3) == 0)
N_imu = N_sens/3;

y = reshape(x, 3, N_imu, []);

end

function g=gravity(lambda,h)
% function g=gravity(lambda,h)
%
% function for calculation of the local gravity vector, in
% the geographic reference frame (same as tangent plane is 
% stationary).
%
% Based upon the WGS_84 Geodetic and Gravity model. For more 
% info see [pp 222-223,1].
%
% lambda -> Latitude [degrees]
% h -> Altitude [m]
% g 
%
% edit: Isaac Skog, 2006-08-17

% degrees to radians
lambda=pi/180*lambda;

gamma=9.780327*(1+0.0053024*sin(lambda)^2-0.0000058*sin(2*lambda)^2);

g=gamma-((3.0877e-6)-(0.004e-6)*sin(lambda)^2)*h+(0.072e-12)*h^2;

g=[0 0 -g]';
return;




function [roll, pitch] = stationaryAcc2rollPitch(u)
%GET_ROLL_PITCH Summary of this function goes here
%   Detailed explanation goes here
f_x=mean(u(1,:));
f_y=mean(u(2,:));
f_z=mean(u(3,:));

roll=atan2(-f_x,f_z);
pitch=atan2(f_y,sqrt(f_x^2+f_z^2));
end

function [roll, pitch] = stationaryAcc2rollPitch_IS(u)
%GET_ROLL_PITCH Summary of this function goes here
%   Detailed explanation goes here
f_u=mean(u(1,:));
f_v=mean(u(2,:));
f_w=mean(u(3,:));

roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));
end

function [y] = triad_mean(x)
%TRIAD_MEAN mean of triad data from 2D matrix
%   y = triad_mean(x)

N_sens = size(x,1);
assert(mod(N_sens,3) == 0)
N_imu = N_sens/3;

y = reshape(mean(reshape(x, 3, N_imu, []),2),3,[]);


end

function [y] = triad_norm(x)
%TRIAD_NORM Summary of this function goes here
%   Detailed explanation goes here

N_sens = size(x,1);
assert(mod(N_sens,3) == 0)
N_imu = N_sens/3;

y = reshape(sqrt(sum(get_triad_form(x).^2, 1)),N_imu,[]);

end

function [w_dot_interp] = w2w_dot_splines(t, w)
%W2W_DOT_SPLINES Interpolate w to w_dot using splines
%
%   w_dot_interp = w2w_dot_splines(t, w)
%
w_dot_interp = zeros(size(w));
for i = 1:3
    pp = spline(t, w(i,:)');
    qq = ppdiff(pp);
    
    w_dot_interp(i,:) = ppval(qq,t);
    
end

end

