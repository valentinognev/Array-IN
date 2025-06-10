function [y] = arrayfun_parfor(f,x)
%ARRAYFUN_PARFOR Summary of this function goes here
%   Detailed explanation goes here
y = zeros(size(x));
tic
parfor i = 1:numel(x)
    y(i) = f(x(i));
end
toc

end

function [X_triu_vec] = get_triu_vec(X)
%GET_TRIU Summary of this function goes here
%   Detailed explanation goes here
X_triu_vec = X(triu(true(size(X))));

end

function [t] = is_any_nan(x)
%IS_ANY_NAN Is any value NaN
t = any(isnan(x),"all");

end

function b = matrix3d2blkdiag(a)

assert(ndims(a) == 3, "Not 3D matrix")

b_cell = num2cell(a, [1,2]);
b = blkdiag(b_cell{:});

endfunction [structB] = mergeStruct(structA,structB)
%mergeStruct Merge two structs
%   If similar then structA has precedence
 f = fieldnames(structA);
 for i = 1:length(f)
    structB.(f{i}) = structA.(f{i});
 end
end

function y = mod1(x, m)
%MOD1		modulo function, but returns m instead of 0
%
% y = mod1(x, m)
%    Return x (mod m), except that if the result is 0, return m instead.
%    This is equal to (x-1 (mod m)) + 1.
%    
%    This function is useful if you have a series of items in a vector v
%    you want to cycle through repeatedly with some index i:
%    use  mod1(i, length(v)).
%
%    Note that mod1, like mod, always returns a positive number.

y = rem(rem(x-1, m) + m, m) + 1;
function show_progress(n,N, N_show)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if mod(n,N_show) == 0
    fprintf("%d / %d\n", n,N);
    waitbar(n/N)
end
 
end

