function [R] = R_bn_spherical(phi,theta, Rb_om, Rn_om)
%R_BN_SPHERICAL from navigation frame to body frame

assert(length(phi) == length(theta), "error");
N = length(phi);
R = zeros(3,3, N);

for n = 1:N
    t = theta(n);
    p = phi(n);
    if nargin<3
        R(:,:,n) = [ cos(t)*sin(p) sin(t)*sin(p)  cos(p);
            cos(t)*cos(p) sin(t)*cos(p) -sin(p);
            -sin(t)        cos(t)         0];
    else
        Rom_n=Rn_om(t)';
        R(:,:,n)=Rb_om(p)*Rom_n;
    end
end


end

