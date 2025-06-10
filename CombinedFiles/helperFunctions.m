function [Y] = commutation_matrix(m,n)
%COMMUTATIONMATRIX Summary of this function goes here
%   Detailed explanation goes here
% [m, n] = size(A);
I = reshape(1:m*n, [m, n]); % initialize a matrix of indices of size(A)
I = I'; % Transpose it
I = I(:); % vectorize the required indices
Y = eye(m*n); % Initialize an identity matrix
Y = Y(I,:); % Re-arrange the rows of the identity matrix
end

function [As, R_square] = compute_As(r_tot)
%COMPUTE_AS The matrix for the rotation 
%
% As = compute_As(r_tot)
% Where r_tot is centered.

assert( all(abs(mean(r_tot,2)) < 10*eps))
K = size(r_tot,2); % Number of acc triads

R_skew = zeros(3,3,K);
for k = 1:K
    R_skew(:,:,k) = skew_sym(r_tot(:,k));
end
R_square = zeros(3,3);
for k = 1:K
    R_square = R_square + R_skew(:,:,k)'*R_skew(:,:,k);
end
As = zeros(3,3,K);
for k = 1:K
    As(:,:,k) = R_square\R_skew(:,:,k);
end

end

function [M] = compute_projection_matrix(H, Q)
%COMPUTE_PROJECTION_MATRIX Summary of this function goes here
%   Detailed explanation goes here
H_t_Q_inv = H'/Q;
M = (H_t_Q_inv * H)\H_t_Q_inv;

end

function [D,L, Q] = factorize_T(T)
%FACTORIZE_T Factorize scale matrix as T = D*L*Q
%   D: Diagonal matrix with scale factors
%   L: Upper triangular matrix with unit diagonal. Account for
%   non-orthogonalities 
%   Q: Rotation matrix 
[R,Q] = rq(T);

S = diag(R);
D = diag(S);
L = D\R;


end

function [m] = get_spherical_motion(t, inp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp(inp.phi, "sinus")
    [phi, phi_dot, phi_dot_2] = get_sinus(t, inp.phi_params);
elseif strcmp(inp.phi, "linear")
    [phi, phi_dot, phi_dot_2] = get_linear(t, inp.phi_params);    
elseif strcmp(inp.phi, "quadratic")
    [phi, phi_dot, phi_dot_2] = get_quadratic(t, inp.phi_params);
elseif strcmp(inp.phi, "poly")    
    [phi, phi_dot, phi_dot_2] = get_polynomial(t, inp.phi_params);
elseif strcmp(inp.phi, "constant")    
    [phi, phi_dot, phi_dot_2] = get_constant(t, inp.phi_params);
else
    error("No correct motion for phi")
end

if strcmp(inp.theta, "sinus")
    [theta, theta_dot, theta_dot_2] = get_sinus(t, inp.theta_params);
elseif strcmp(inp.theta, "linear")
    [theta, theta_dot, theta_dot_2] = get_linear(t, inp.theta_params);
elseif strcmp(inp.theta, "quadratic")
    [theta, theta_dot, theta_dot_2] = get_quadratic(t, inp.theta_params);
elseif strcmp(inp.theta, "poly")
    [theta, theta_dot, theta_dot_2] = get_polynomial(t, inp.theta_params);
elseif strcmp(inp.theta, "constant")
    [theta, theta_dot, theta_dot_2] = get_constant(t, inp.theta_params);
else
    error("No correct motion for theta")
end
m.phi = phi;
m.phi_dot = phi_dot;
m.phi_dot_2 = phi_dot_2;

m.theta = theta;
m.theta_dot = theta_dot;
m.theta_dot_2 = theta_dot_2;

end

function [s, s_dot, s_dot_2] = get_constant(t, inp)
    if isfield(inp,"A")
        A = inp.A;
    else
        A = 1;
    end
    s = A*ones(size(t));
    s_dot = zeros(size(t));
    s_dot_2 = zeros(size(t));

end
function [s, s_dot, s_dot_2] = get_sinus(t, inp)
    if isfield(inp,"A") 
        A = inp.A;
    else
        A = 1;
    end
    if isfield(inp,"f") 
        f = inp.f;
    else
        f = 1;
    end
    if isfield(inp,"b") 
        b = inp.b;
    else
        b = 0;
    end
    s = A.*sin(2*pi*f*t) + b;
    s_dot = A.*cos(2*pi*f*t)*2*pi*f;
    s_dot_2 = -A.*sin(2*pi*f*t)*(2*pi*f)^2;    
end
function [s, s_dot, s_dot_2] = get_linear(t, inp)
    if isfield(inp,"A") 
        A = inp.A;
    else
        A = 1;
    end
    s = A*t;
    s_dot = A*ones(size(t));
    s_dot_2 = zeros(size(t));
end
function [s, s_dot, s_dot_2] = get_quadratic(t, inp)
    if isfield(inp,"A") 
        A = inp.A;
    else
        A = 1;
    end
    s = A*t.^2;
    s_dot = 2*A*t;
    s_dot_2 = 2*A*ones(size(t));
end
function [s, s_dot, s_dot_2] = get_polynomial(t, inp)
    p = inp.p;
    s = polyval(p, t);
    p1 = polyder(p);
    s_dot = polyval(p1, t);
    p2 = polyder(p1);
    s_dot_2 = polyval(p2, t);
end
function [s,l,q] = get_T_components(T)
%GET_T_COMPONENTS Get the components of the T matrix
%   s: scale factors 
%   l: angles for non-orthogonalities
%   q: rotation vector 
[S,L,Q] = factorize_T(T);

s = diag(S);

l = zeros(3,1);
l(1) = L(1,2);
l(2) = L(1,3);
l(3) = L(2,3);

q = logSO3(Q);
end

function [R] = initialAttitude2Rotm(u)
%INITIAL_ATTITUDE Summary of this function goes here
%   Detailed explanation goes here

f_u=mean(u(1,:));
f_v=mean(u(2,:));
f_w=mean(u(3,:));

roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));
heading = 0;

R = eul2rotm([roll pitch heading], 'XYZ');

end

function [norm_v] = norm_time(x)
%NORM_TIME Summary of this function goes here
%   Detailed explanation goes here

norm_v = sqrt(sum(x.^2,1));
end

function qq = ppdiff(pp,j)
%PPDIFF Differentiate piecewise polynomial.
%   QQ = PPDIFF(PP,J) returns the J:th derivative of a piecewise
%   polynomial PP. PP must be on the form evaluated by PPVAL. QQ is a
%   piecewise polynomial on the same form. Default value for J is 1.
%
%   Example:
%       x = linspace(-pi,pi,9);
%       y = sin(x);
%       pp = spline(x,y);
%       qq = ppdiff(pp);
%       xx = linspace(-pi,pi,201);
%       plot(xx,cos(xx),'b',xx,ppval(qq,xx),'r')
%
%   See also PPVAL, SPLINE, SPLINEFIT, PPINT

%   Author: Jonas Lundgren <splinefit@gmail.com> 2009

if nargin < 1, help ppdiff, return, end
if nargin < 2, j = 1; end

% Check diff order
if ~isreal(j) || mod(j,1) || j < 0
    msgid = 'PPDIFF:DiffOrder';
    message = 'Order of derivative must be a non-negative integer!';
    error(msgid,message)
end

% Get coefficients
coefs = pp.coefs;
[m, n] = size(coefs);

if j == 0
    % Do nothing
elseif j < n
    % Derivative of order J
    D = [n-j:-1:1; ones(j-1,n-j)];
    D = cumsum(D,1);
    D = prod(D,1);
    coefs = coefs(:,1:n-j);
    for k = 1:n-j
        coefs(:,k) = D(k)*coefs(:,k);
    end
else
    % Derivative kills PP
    coefs = zeros(m,1);
end

% Set output
qq = pp;
qq.coefs = coefs;
qq.order = size(coefs,2);
function [dtheta] = q_minus(q1,q2)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
% Half the angle 
o1 = ones(length(q1),1);
o1(parts(q1) < 0) = -1;
q1 = q1.*o1;

o2 = ones(length(q2),1);
o2(parts(q2) < 0) = -1;
q2 = q2.*o2;

dtheta = compact(log(conj(q1).*q2)).*2;
assert(all(abs(dtheta(:,1)) < 100*eps));
dtheta = dtheta(:,2:4); %

end

function [R, Q] = rq(T)
%RQ RQ factorization
%   Same as QR and R have positive diagoanls 

[Q,~] = qr(flipud(T)');
Q = fliplr(Q); % Upper triangularize T{1} from the left
Q = Q*diag(diag(sign(Q))); % To not change coordinate system orientation
Q = Q';

R = T*Q';

end

function A=skew_sym(a)

A=[   0  -a(3)  a(2); 
    a(3)    0  -a(1); 
   -a(2)  a(1)    0];


endfunction [R] = solve_Wahbas_problem(W,V)
%solve_Wahbas_problem Estimate initial rotation matrix from gravity 
%   R = argmin sum_{k,n} || w_{k,n} - R*v_{k,n}||^2
%  where R in SO(3).
assert(all(size(W) == size(V)))
N = size(W,2);
K = size(W,1)/3;
B = zeros(3);
inds = reshape(1:3*K,3,[]);
for n = 1:N
    for k = 1:K
        kk = inds(:,k);
        B = B + W(kk,n)*V(kk,n)';
    end
end

[U,~,V] = svd(B);
M = eye(3);
M(3,3) = det(U)*det(V);
R = U*M*V';

end

function [D,V,P_norm] = stochastic_observability(P)
%STOCHASTIC_OBSERVABILITY Summary of this function goes here
%   Detailed explanation goes here

assert(ndims(P) == 3)

P0 = P(:,:,1);
n = size(P,1);

assert(isdiag(P0))

F = inv(sqrt(P0));

P_norm = zeros(size(P));


for k = 1:size(P,3)
    P_k_in = P(:,:,k);
    
    P_k_in_1 = F*P_k_in*F;

    % Normalize to unit norm for the eigen values 
    P_norm_k = P_k_in_1./trace(P_k_in_1);
    
    P_norm(:,:,k) = P_norm_k;
    
    
end

% eigenshuffle: Consistent sorting for an eigenvalue/vector sequence
% [Vseq,Dseq] = eigenshuffle(Asequence)

[V, D] = eigenshuffle(P_norm);


end

