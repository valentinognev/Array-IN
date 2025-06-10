function [out] = myfminunc(objfun, x0, options)

[x_opt,fval,exitflag,output,grad,hessian] = fminunc(objfun,x0,options);

out = struct;
out.x_opt = x_opt;
out.fval = fval;
out.exitflag = exitflag;
out.output = output;
out.grad = grad;
out.hessian = hessian;



endfunction [h, g] = numerical_hessian_forward_diff(f,x)

epsilon = 1e-5; 
epsilon_inv = 1/epsilon;

nx = length(x); % Dimension of the input x;
f0 = feval(f, x); % caclulate f0, when no perturbation happens

f_e = zeros(nx,1);
g = zeros(nx,1);
% Do perturbation
for i = 1:nx
    x_ = x;
    x_(i) =  x(i) + epsilon;
    f_e(i) = feval(f, x_);
    
    g(i) = (f_e(i) - f0) .* epsilon_inv;
end

Boolind = triu(true(nx,nx));
Boolind = Boolind(:);

n_2e = nx*(nx + 1)/2;
f_2e_vec = zeros(n_2e, 1);

inds = triu(reshape(1:nx^2,nx,nx));
inds_vec = inds(Boolind);

for n = 1:n_2e
    ind = inds_vec(n);
    [i,j] = ind2sub([nx,nx], ind);

    x_ = x;
    if i == j
        x_(i) =  x(i) + 2*epsilon;
    else
        x_(i) =  x(i) + epsilon;
        x_(j) =  x(j) + epsilon;
    end
    f_2e_vec(n) = feval(f, x_);
end

f_2e = zeros(nx, nx);
f_2e(Boolind) = f_2e_vec;

f_2e = f_2e + triu(f_2e,1)';
h = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        h(i,j) = (f_2e(i,j) - f_e(i) - f_e(j)  + f0)*epsilon_inv^2;
    end
end

endfunction jac = numeric_jacobian(f, x)
% Calculate Jacobian of function f at given x
epsilon = 1e-6; 
epsilon_inv = 1/epsilon;
nx = length(x); % Dimension of the input x;
f0 = feval(f, x); % caclulate f0, when no perturbation happens
nf = length(f0);
jac = zeros(nf,nx);
% Do perturbation
for i = 1:nx
    x_ = x;
    x_(i) =  x(i) + epsilon;
    jac(:, i) = (feval(f, x_) - f0) .* epsilon_inv;
endfunction [h, g] = pnumerical_hessian_forward_diff(f,x)

epsilon = 1e-5; 
epsilon_inv = 1/epsilon;

nx = length(x); % Dimension of the input x;
f0 = feval(f, x); % caclulate f0, when no perturbation happens

f_e = zeros(nx,1);
g = zeros(nx,1);
% Do perturbation
parfor i = 1:nx
    x_ = x;
    x_(i) =  x(i) + epsilon;
    f_e(i) = feval(f, x_);
    
    g(i) = (f_e(i) - f0) .* epsilon_inv;
end

Boolind = triu(true(nx,nx));
Boolind = Boolind(:);

n_2e = nx*(nx + 1)/2;
f_2e_vec = zeros(n_2e, 1);

inds = triu(reshape(1:nx^2,nx,nx));
inds_vec = inds(Boolind);

parfor n = 1:n_2e
    ind = inds_vec(n);
    [i,j] = ind2sub([nx,nx], ind);

    x_ = x;
    if i == j
        x_(i) =  x(i) + 2*epsilon;
    else
        x_(i) =  x(i) + epsilon;
        x_(j) =  x(j) + epsilon;
    end
    f_2e_vec(n) = feval(f, x_);
end

f_2e = zeros(nx, nx);
f_2e(Boolind) = f_2e_vec;

f_2e = f_2e + triu(f_2e,1)';
h = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        h(i,j) = (f_2e(i,j) - f_e(i) - f_e(j)  + f0)*epsilon_inv^2;
    end
end

endfunction [out] = runfmincon(objfun, x0, options, varargin)

opt = struct;
if ~isempty(varargin) 
    assert(mod(length(varargin),2) == 0)
    names = varargin{1:2:end};
    values = varargin{2:2:end};
    for k = 1:length(names)
        opt.(names{k}) = values{k};
    end
end
% Set up shared variables with outfun
history.x = [];
history.fval = [];
searchdir = [];
x0 = reshape(x0,[],1);
% Call optimization
% x0 = [-1 1];
% options = optimoptions(@fmincon,'OutputFcn',@outfun,...
%     'Display','iter','Algorithm','active-set');
options.OutputFcn = @outfun;
if isfield(opt,"hessian_and_grad") && opt.("hessian_and_grad")
    printf("Save final hessian and gradient\n")
    [x_opt,fval,exitflag,output,grad,hessian] = fminunc(objfun,x0,options);
else
    [x_opt,fval,exitflag,output] = fminunc(objfun,x0,options);
end
history.x = reshape(history.x, length(x0), []);
out = struct;
out.x_opt = x_opt;
out.fval = fval;
out.exitflag = exitflag;
out.output = output;
if isfield(opt,"hessian_and_grad") && opt.("hessian_and_grad")
    out.grad = grad;
    out.hessian = hessian;
end
out.history = history;
out.searchdir = searchdir;



    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                % hold on
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x; x];
                % Concatenate current search direction with
                % searchdir.
                searchdir = [searchdir;...
                    optimValues.searchdirection'];
                %  plot(x(1),x(2),'o');
                % Label points with iteration number and add title.
                % Add .15 to x(1) to separate label from plotted 'o'.
                % text(x(1)+.15,x(2),...
                %     num2str(optimValues.iteration));
                % title('Sequence of Points Computed by fmincon');
            case 'done'
                % hold off
            otherwise
        end
    end
endfunction hf =solveHessian(test_function, a)
% Objective: Generates Hessian of a function at some point
%-----------------------------------------------------------------------
% f=solveHessian(a,test_function)
% where a=input vector
%       test_function=objective function
%-----------------------------------------------------------------------
% Output: f= Hesian matrix
%-----------------------------------------------------------------------

% Code by:
% Salil Sharma
% May 3, 2017
%-----------------------------------------------------------------------

l=length(a); %Hessian would be lxl matrix
ep=0.0001; % step size for numerical diffrentiation
valf=test_function(a); % value of obj function at a
ep2=ep*ep;
ep3=4*ep*ep;
hf = zeros(l,l);
for i=1:length(a)
    x1=a;
    x1(i)=a(i)-ep; %Change ith element in x1
    x2=a;
    x2(i)=a(i)+ep; %Change ith element in x2
    hf(i,i)=(test_function(x2)-2*valf+test_function(x1))/ep2; % diagonal entries
    j=i+1;
    while j<=length(a) % Loop computes the rest of the elements of the Hessian matrix
        x1(j)=a(j)-ep; % Lower the value of step size
        x2(j)=a(j)+ep; % Increment the value of step size
        v4=test_function(x1); % compute the respective values
        v1=test_function(x2); % compute the respective values
        x1(j)=x1(j)+2*ep; 
        x2(j)=x2(j)-2*ep;
        v2=test_function(x1);
        v3=test_function(x2);
        hf(i,j)=(v1+v4-v2-v3)/ep3;
        hf(j,i)=hf(i,j); % d2f/dxdy is same as that of d2f/dydx
        x1(j)=a(j);
        x2(j)=a(j); 
        j=j+1;
    end
end

% for i=1:length(a)    
%     while j<=length(a) % Loop computes the rest of the elements of the Hessian matrix
%         hf(j,i)=hf(i,j); % d2f/dxdy is same as that of d2f/dydx
%     end
% end
end
function d = spline_derivative(t, w, j)

if nargin < 3
    j = 1;
end
d = zeros(size(w));
for i = 1:3
    pp = spline(t, w(i,:)');
    qq = ppdiff(pp,j);    
    d(i,:) = ppval(qq,t);    
end

end

