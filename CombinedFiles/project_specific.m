function [err] = calculate_trajectory_error(S_hat,S_true,varargin)
%CALCULATE_TRAJECTORY_ERROR Summary of this function goes here
%   Detailed explanation goes here
if nargin > 2
    section = varargin{1};
else
    section = 1:size(S_true.R,3);
end
err = struct;
err.R = errorSO3(S_hat.R, S_true.R(:,:,section));
err.R_deg = rad2deg(err.R);
if isfield(S_hat, "w")
    err.w = S_hat.w - S_true.w(:,section); 
    err.w_deg = rad2deg(err.w);
end
err.p = S_hat.p - S_true.p(:,section);
err.v = S_hat.v - S_true.v(:,section);

if isfield(S_hat, "b_a")
    inds_bias = 1:size(S_hat.b_a,1);
    err.b_a = S_hat.b_a - S_true.b_a(inds_bias,section);    
end
if isfield(S_hat, "b_g")
    err.b_g = S_hat.b_g - S_true.b_g(:,section); 
    err.b_g_deg = rad2deg(err.b_g);
end

end

function [Q_u] = compensate_covariance(Q_y,T)
%COMPENSATE_MEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

L = chol(Q_y, "lower");

q = T\L;

Q_u = q*q';

end

function [u] = compensate_measurements(y,T,b)
%COMPENSATE_MEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
if ndims(T) == 3
    T_diag = matrix3d2blkdiag(T);
else
    T_diag = T;
end
b = reshape(b,[],1);
u = T_diag\(y - b);
end

function [err] = compute_error(S,S_ref)
%COMPUTE_ERROR Summary of this function goes here
%   Detailed explanation goes here
err = struct;
if isfield(S_ref,"R") && isfield(S,"R")
    try
        err.R = errorSO3(S.R, S_ref.R);
    catch
        warning('Angle Error is too high.');
    end
end

if isfield(S_ref,"v") && isfield(S,"v")
    err.v = S.v - S_ref.v;
end

if isfield(S_ref,"p") && isfield(S,"p")
    err.p = S.p - S_ref.p;
end

if isfield(S_ref,"w") && isfield(S,"w")
    err.w = S.w - S_ref.w;
end

if isfield(S_ref,"omega_dot") && isfield(S,"omega_dot")
    err.omega_dot = S.omega_dot - S_ref.omega_dot;
end

if isfield(S_ref,"v_dot") && isfield(S,"v_dot")
    err.v_dot = S.v_dot - S_ref.v_dot;
end

if isfield(S_ref,"s") && isfield(S,"s")
    err.s = S.s - S_ref.s;
end

if isfield(S_ref,"b_g") && isfield(S,"b_g")
    err.b_g = S.b_g - S_ref.b_g;
end

if isfield(S_ref,"b_s") && isfield(S,"b_s")
    err.b_s = S.b_s - S_ref.b_s;
end

if isfield(S_ref,"T_a") && isfield(S,"T_a")
    err.T_a = S.T_a - S_ref.T_a;
end

if isfield(S_ref,"b_omega_dot") && isfield(S,"b_omega_dot")
    err.b_omega_dot = S.b_omega_dot - S_ref.b_omega_dot;
end

end

function [T,b] = estimate_T_and_b(y,u,Q)
%ESTIMATE_T_AND_B Summary of this function goes here
%   Detailed explanation goes here
assert(all(size(y) == size(u)))
A = zeros(12,12);
b = zeros(12,2);
for n = 1:size(y,2)
    H_n = [kron(u(:,n)',eye(3)) eye(3)];
    Ht_Q_inv = H_n'/Q;
    A = A + Ht_Q_inv*H_n;
    b = b + Ht_Q_inv*y(:,n);
end
Tb = A\b;

T = reshape(Tb(1:9),3,3);
b = Tb(10:12);

end

function res = get_release_inds(time, release_times, IN_time, T)


%IN_time = 5; 
IN_samples = IN_time/T;
inds_growth = zeros(IN_samples, length(release_times));
for i_t = 1:length(release_times)
    inds_t = find(time >= release_times(i_t));
    inds_growth(:,i_t) = inds_t(1:IN_samples);
    
end
IN_time_array = (0:IN_samples-1)*T;
res = struct;
res.inds_growth = inds_growth;
res.IN_time_array = IN_time_array;

endfunction [S] = interpolate_pos_and_rotation(S,imu_time, rig_time)
%INTERPOLATE_POS_AND_ROTATION Interpolate IMU pos and rotation estimates 
% to rig time
S.p_rig_time = zeros(3, length(rig_time));
for i = 1:3
    S.p_rig_time(i,:) = interp1(imu_time, S.p(i,:), rig_time, "pchip");
end


% Find the fractional indices using linear interpolation 
inds_time = interp1(imu_time, 1:length(imu_time), rig_time, "linear");
R_rig_time = zeros(3,3,length(inds_time));
for n = 1:length(inds_time)
    
    frac = inds_time(n) - floor(inds_time(n));
    if frac > 1
        warning("fraction larger than 1")
    elseif frac == 0
        % Same point in time
        R_rig_time(:,:,n) = S.R(:,:, round(inds_time(n)));
    else
        % Calculate the rotation vector and scale it
        left_ind = floor(inds_time(n));
        right_ind = left_ind + 1;
        R_left = S.R(:,:, left_ind);
        R_right = S.R(:,:, right_ind);
        
        theta = logSO3(invSO3(R_left)*R_right);
        R_rig_time(:,:,n) = R_left*expSO3(frac*theta);
    end
end
S.R_rig_time = R_rig_time;

end

function [u, Qu] = lsq_triad(y,Q)
%LSQ_TRIAD Weighted Mean of triad
%   y = Hu + e , e ~ N(0,Q)
%  u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
% Where u is triad
assert(mod(size(y,1),3) == 0)
assert(size(Q,1) == size(Q,2))

N = size(y,1)/3;
H = repmat(eye(3), N, 1);
L = chol(Q, "lower");

t1 = H'/L;
t2 = (L')\y;
u = (t1 * t1')\(t1*t2);
Qu = inv(t1 * t1');

end

function [u, Qu] = lsq_triad_naive(y,Q)
%LSQ_TRIAD Summary of this function goes here
%   y = Hu + e , e ~ N(0,Q)
%  u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
% Where u is triad
assert(mod(size(y,1),3) == 0)
assert(size(Q,1) == size(Q,2))

N = size(y,1) / 3;
H = kron(ones(N,1),eye(3));
HT_Q_inv = (H')/Q;

u = (HT_Q_inv * H)\(HT_Q_inv * y);

Qu = inv(HT_Q_inv * H);

end



function [S_out] = rotate_measurements(S_in,R)
%ROTATE_MEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
S_out = struct;
S_out.y = R*S_in.y;
S_out.Q = R*S_in.Q*R';
S_out.Q_inv = inv(S_out.Q);
end

function res = run_filter(sensorData, initData, my_settings, S_ref, myFilter)
%RUN_FILTER Run filter and calculate error 
%   Detailed explanation goes here
res = struct;
[res.filt, res.pred] = myFilter(sensorData, initData, my_settings);
err = struct;
if isfield(S_ref,"R")
    try
        err.R = errorSO3(res.filt.R, S_ref.R);
    catch
        warning('Angle Error is too high.');
    end

end

if isfield(S_ref,"v")
    err.v = res.filt.v - S_ref.v;
end

if isfield(S_ref,"p")
    err.p = res.filt.p - S_ref.p;
end

if isfield(S_ref,"w") && isfield(res.filt,"w")
    err.w = res.filt.w - S_ref.w;
end

if isfield(S_ref,"omega_dot")
    err.omega_dot = res.pred.omega_dot - S_ref.omega_dot;
end

if isfield(S_ref,"v_dot")
    err.v_dot = res.pred.v_dot - S_ref.v_dot;
end

res.err = err;

end


function res = run_filter_w_error(sensorData, initData, my_settings, S_ref, myFilter)
%RUN_FILTER Run filter and calculate error 
%   Detailed explanation goes here
res = struct;
[res.filt, res.pred] = myFilter(sensorData, initData, my_settings);
err = struct;
err.filt = compute_error(res.filt, S_ref);
err.pred = compute_error(res.pred, S_ref);
res.err = err;

end


