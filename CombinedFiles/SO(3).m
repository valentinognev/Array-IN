function adw = adjSO3(w)

adw = HatSO3(w);function AdR = AdSO3(R)

AdR = R;
function [ x_dual_hat ] = dualHatSO3( x )
%w_hat*x = x_dual_hat*w

x_dual_hat = - HatSO3(x);


end

function [err] = errorSO3(Rhat, Rtrue)
%UNTITLED Summary of this function goes here
%   Need to be same length 

N = size(Rhat,3);

%% Calculate errors
err = zeros(3, N);
for n = 1:N
    if ndims(Rtrue) == 3
        err(:, n) = logSO3(invSO3(Rhat(:,:,n))*Rtrue(:,:,n));
    else
        err(:, n) = logSO3(invSO3(Rhat(:,:,n))*Rtrue);
    end
end

end

function R = expSO3(w)

normw = norm(w);

if(normw == 0)
    R = eye(3);
    return;
end
w_hat = HatSO3(w);

R = eye(3) + sin(normw)*w_hat/normw + (1-cos(normw))*w_hat*w_hat/(normw^2);

function w_hat = HatSO3(w)

w_hat = [0 -w(3) w(2);...
         w(3) 0 -w(1);...
         -w(2) w(1) 0];
     
% w_hat = zeros(3,3);
% w_hat(2,1) = w(3);
% w_hat(3,1) = -w(2);
% 
% w_hat(1,2) = -w(3);
% w_hat(3,2) = w(1);
% 
% w_hat(1,3) = w(2);
% w_hat(2,3) = -w(1);function [Rinv,errorFlag] = invSO3(R)

errorFlag = 0;
Rinv = R';
function [w,errorflag] = logSO3(R)

phy = acos((trace(R)-1)/2);
if(abs(phy)> pi)
    error('angle sup�rieur � pi');
end



if(phy == 0)
    w = zeros(3,1);
elseif (abs(phy) == pi)
    
    A = (R-eye(3))/2;
    w1 = sqrt(-((A(2,2) + A(3,3) - A(1,1))/2));
    
    
    w2 = sqrt(-((A(1,1) + A(3,3) - A(2,2))/2));
    w3 = sqrt(-((A(1,1) + A(2,2) - A(3,3))/2));
    
    if(w1~=0)
        
        if(A(1,2) < 0)
            w2 = -w2;
        end
        if(A(1,3) < 0)
            w3 = -w3;
        end
        
    elseif(w2~=0)
        
        if(A(2,3) < 0)
            w3 = -w3;
        end
    end
    
    w = [w1;w2;w3]*phy;
    
else
    w_hat = (R-R.')/(2*sin(phy))*phy;%on remultiplie par phy pour retrouver le vecteur avec sa norme originale
    w = VecSO3(w_hat);
end

errorflag = 0;
endfunction Rnorm = normalizeSO3(R)

[u,s,v] = svd(R);

Rnorm = u*v';
endfunction Phiw = PhiSO3(w)
% Left-Jacobian to SO(3)
% sum_k 1/(k + 1)! ad(w)^k
normw = norm(w);

% if(normw > pi/2)
%     error('formula not sure')
% end
if(normw > 0)
    adw = adjSO3(w);
    
    Phiw = eye(3) + (1/(2*normw^2))*(4-normw*sin(normw)-4*cos(normw))*adw+...
        (1/(2*normw^3))*(4*normw-5*sin(normw)+normw*cos(normw))*adw^2+...
        (1/(2*normw^4))*(2-normw*sin(normw)-2*cos(normw))*adw^3+...
        (1/(2*normw^5))*(2*normw-3*sin(normw)+normw*cos(normw))*adw^4;
else
    Phiw = eye(3);
end
end
function w = VecSO3(w_hat)

w = [w_hat(3,2); w_hat(1,3); w_hat(2,1)];
