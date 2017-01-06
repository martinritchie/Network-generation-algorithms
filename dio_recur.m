function [solution, rl] = dio_recur(coeffs, sol,rl)
%% Description (I)
% This function recursively solves diophantine equations in n terms. Ensure
% that global variables are correctly cleared between successive calls of
% this function.
%% Input / Output
%   coeffs: coefficients of the diophantine equation (enter in ascending
%           order),
%      sol: the constant term,
%       rl: recursion level,
% solution: solution space corresponding to the constant term
%% Example call
% solution = dio_recur([1 2 2 3], 5)
%% Procedure
set(0,'RecursionLimit',500)
global sol_space index
% The following initialises a few key variables on the first call.
if nargin < 3
    rl = 0;
    % index: index of the current solution within the solution space.
    index = 1;
else
    rl = rl + 1;
end

if length(coeffs) ==1
    sol_space = round(sol/coeffs);
elseif length(coeffs) > 2
    dim1 = length(coeffs);
    q = floor(sol./coeffs(end));
    for l = 0:q
        % when comming out of a recursion we duplicate the preceding parent
        % lines of the solution space.
        if l > 0
            solution_space_rep(index);
        end
        solution_space(l,index,dim1);
        dio_recur(coeffs(1:end-1), sol-l*coeffs(end),rl);
    end
else
    q = floor(sol./coeffs);
    count2 = 0;
    for j = 0:q(2)
        for i = 0:q(1)
            if i*coeffs(1) + j*coeffs(2) ==sol
                % when comming out of a recursion we duplicate the
                % preceding parent lines of the solution space.
                if count2 >  0
                    solution_space_rep(index);
                end
                count2 = count2 + 1;
                solution_space(i,index, 1);
                solution_space(j,index, 2);
                index = index + 1;
            end
        end
    end
end
solution = sol_space;
end
function solution_space(x,index,index2)
%% Description (II)
% This function keeps track of, and update the solution space.
global sol_space k
sol_space(index, index2)  = x;
k = k + 1;
end
function solution_space_rep(index)
%% Description (III)
% This function is used to repeat the previous part of the solution space.
global sol_space
sol_space(index, 3:end)  = sol_space(index-1, 3:end);
end