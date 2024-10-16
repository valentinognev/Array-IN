function [A] = compute_A_non_center(sensorPos, T)
%COMPUTE_AS The matrix for the rotation 
% Based on: Appendix A of "Inertial Navigation using an Inertial sensor array"
%
% As = compute_As(sensorPos)
% Where sensorPos is centered - position relative to centroid.

K = size(sensorPos,2); % Number of acc triads

R_skew = cell(K,1);
for k = 1:K
    R_skew{k} = skew_sym(sensorPos(:,k));
end
H = [-cat(1,R_skew{:}), repmat(eye(3),K,1)];

if nargin == 2
    H = matrix3d2blkdiag(T)*H;
end

A = (H'*H)\H';  % 5b - Inertial Navigation using an Inertial sensor array



