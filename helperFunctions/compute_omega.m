function vels = compute_omega(inp)
theta = inp.theta;
theta_dot = inp.theta_dot;
theta_dot_2 = inp.theta_dot_2;

phi = inp.phi;
phi_dot = inp.phi_dot;
phi_dot_2 = inp.phi_dot_2;

N = length(theta);

omega=[phi_dot*0;
       phi_dot;
       theta_dot];
omegadot=[-phi_dot.*theta_dot;
           phi_dot_2
           theta_dot_2];
Rn_om = @(theta) [cos(theta), -sin(theta), 0; ...
                  sin(theta),  cos(theta), 0;
                       0    ,       0    , 1];
Rb_om = @(phi)   [sin(phi)  ,      0     ,  cos(phi);
                  cos(phi)  ,      0     , -sin(phi);
                     0      ,      1     ,     0    ];

w_n = zeros(3,N);
w_dot_n = zeros(3,N);
for i=1:length(theta)
    w_n(:,i)=Rn_om(theta(i))*omega(:,i);
    w_dot_n(:,i) = Rn_om(theta(i))*omegadot(:,i);
end
% % The angular velocity in navigation frame
% w_n = zeros(3,N);
% w_n(1,:) = phi_dot.*(-sin(theta)); % x
% w_n(2,:) = phi_dot.*(cos(theta));  % y
% w_n(3,:) = theta_dot;  % z
% 
% % The angular acceleration in navigation frame
% w_dot_n = zeros(3,N);
% w_dot_n(1,:) = -phi_dot_2.*sin(theta) - phi_dot.*theta_dot.*cos(theta); % x
% w_dot_n(2,:) = phi_dot_2.*cos(theta) - phi_dot.*theta_dot.*sin(theta);  % y
% w_dot_n(3,:) = theta_dot_2;  % z

% Body frame
w_b = zeros(3,N);
w_dot_b = zeros(3,N);
for i=1:length(phi)
    w_b(:,i) = Rb_om(phi(i))*omega(:,i);
    w_dot_b(:,i) = Rb_om(phi(i))*omegadot(:,i);
end
% w_b = zeros(3,N);
% w_b(1,:) =  theta_dot.*cos(phi); % r
% w_b(2,:) = -theta_dot.*sin(phi); % phi
% w_b(3,:) =  phi_dot; % theta
% 
% w_dot_b = zeros(3,N);
% w_dot_b(1,:) =  theta_dot_2.*cos(phi) - phi_dot.*theta_dot.*sin(phi); % r
% w_dot_b(2,:) = -theta_dot_2.*sin(phi) - phi_dot.*theta_dot.*cos(phi); % phi
% w_dot_b(3,:) =  phi_dot_2; % theta

vels.w_b = w_b;
vels.w_dot_b = w_dot_b;

vels.w_n = w_n;
vels.w_dot_n = w_dot_n;
end
