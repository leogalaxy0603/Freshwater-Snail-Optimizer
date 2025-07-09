function [bestSolution, bestFitness] = snail_optimize(func, dim, lb, ub, maxIter, popSize)
% Freshwater Snail Optimizer
% This function implements single-objective continuous optimization based on the behavior pattern of freshwater snails moving in water flow.
% Calling format:
% [bestSolution, bestFitness] = snail_optimize(func, dim, lb, ub, maxIter, popSize)
% Parameters:
% func - objective function handle (@function name) or function name string. The objective function should accept a 1×dim vector input and return a scalar fitness value (the smaller the better).
% dim - the dimension of the problem (the number of decision variables).
% lb - the lower bound of the decision variable (scalar or 1×dim vector).
% ub - the upper bound of the decision variable (scalar or 1×dim vector).
% maxIter - the maximum number of iterations.
% popSize - population size (the number of individuals in the snail population).
% Output:
% bestSolution - the best solution obtained by the search (1×dim vector).
% bestFitness - the objective function value corresponding to the best solution.
% If you have any questions, please contact us
%
%         Email：gaoyuanliu@sanyau.edu.cn
%
%----------The source code will be released uponacceptance--------------