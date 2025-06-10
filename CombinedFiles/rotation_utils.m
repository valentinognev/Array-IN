function [R_mean, list_r] = average_rotation(Rs, nb_it_max, tol_r)
%average_rotation Average rotation from multiple rotation matrices
%   [R_mean, list_r] = average_rotation(Rs, nb_it_max, tol_r)
%   size(Rs) == [3,3,N], N number of rotations
%   nb_it_max: max number of iterations, default 20 
%   tol_r: tolarance for deviations, default 1e-10
%   R_mean: Average rotation 
%   list_r: residuals in lie algebra 

if nargin == 1
    nb_it_max = 20;
    tol_r     = 1e-10;   % [1]
end
number_rotations = size(Rs,3);
R_mean = Rs(:,:,1); % First approx of R [1]
for nb_it = 1:nb_it_max % [2]
    list_r = nan(3,number_rotations);  % [3]
    for i = 1:number_rotations
        list_r(:,i) = logSO3(R_mean'*Rs(:,:,i));
    end
    r = mean(list_r,2);
    
    fprintf("%d/%d: tol: %.3e / %.3e\n", nb_it, nb_it_max, norm(r), tol_r);
    if norm(r) < tol_r % [4]
        break
    end
    R_mean = R_mean * expSO3(r); % Update [7]
    
end % [8]
if nb_it == nb_it_max
    error('the maximum number of iteration where reached')
end


end

function R_inter = interpolate_rotation(t, R, t_inter)
%INTERPOLATE_POS_AND_ROTATION Interpolate IMU pos and rotation estimates 
% R_inter = interpolate_rotation(t, R, t_inter)
% t and R are data points 
% t_inter is the time points where interpolation should occur
% R_inter is the interpolated rotation matrix

S = size(R);
assert(length(t) == S(3))
assert(all(S(1:2) == [3 3]))


% Find the fractional indices using linear interpolation 
inds_imu_time = interp1(t, 1:length(t), t_inter, "linear");
R_inter = NaN(3,3,length(inds_imu_time));
for n = 1:length(inds_imu_time)
    % Extrapolation set to NaN
    if isnan(inds_imu_time(n))
        continue
    end
    
    frac = inds_imu_time(n) - floor(inds_imu_time(n));
    if frac > 1
        warning("fraction larger than 1")
    elseif frac == 0
        % Same point in time
        R_inter(:,:,n) = R(:,:, round(inds_imu_time(n)));
    else
        % Calculate the rotation vector and scale it
        left_ind = floor(inds_imu_time(n));
        right_ind = left_ind + 1;
        R_left = R(:, :, left_ind);
        R_right = R(:, :, right_ind);
        
        theta = logSO3(invSO3(R_left)*R_right);
        R_inter(:,:,n) = R_left*expSO3(frac*theta);
    end
end


end

function [E] = my_rotm2eul(R)
%MY_ROTM2EUL Rotation matrix to euler angles [roll, pitch, yaw]
%
%   E = my_rotm2eul(R)
%
%   Roll: around x-axis
%   Pitch: around y-axis
%   Yaw: around z-axis (heading)
%   R (3,3,N) ->  E (3, N) 

% rotm2eul gives [yaw, pitch, roll] intrinsic rotation
% R = R_z(yaw)*R_y(pitch)*R_z(roll)
% unwrap: adds 2pi when wrapping 
% flipud to get in order [roll, pitch, yaw]
E = flipud(unwrap(rotm2eul(R, 'ZYX'))'); 

end

function w = R2w_central_diff(R,t)
%R2W_CENTRAL_DIFF Rotation matrix 2 angular velocity using central 
%difference 
% 
%   w = R2w_central_diff(R,t)
%
%   Based on:
%   R_{t+1} = R_{t} exp_SO3(w*t)
%   w in body frame

w = nan(3,length(t));

for n = 2:length(t) - 1
    dt = t(n+1) - t(n-1);
    w(:,n) = logSO3(R(:,:,n-1)'*R(:,:,n+1))/dt;
end


end

function [R] = rotationMatrixFromTwoUnitVectors(a,b)
%rotationMatrixFromTwoUnitVectors Find rotation matrix from a to b
%   Detailed explanation goes here
a = a./norm(a);
b = b./norm(b);

v = cross(a,b);
s = norm(v);
c = dot(a,b);

R = eye(3) + skew_sym(v) + skew_sym(v)^2 *(1-c)/s^2;

end

