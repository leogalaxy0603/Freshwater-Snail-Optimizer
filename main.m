clear all
clc
close all

dim = 50;                % Problem dimensionality: number of decision variables
lb = -100;               % Lower bound: scalar means every dimension has lower bound -100
ub = 100;                % Upper bound: scalar means every dimension has upper bound 100
maxIter = 1000;           % Maximum number of iterations/generations for the optimizer
popSize = 50;            % Population size: number of candidate solutions

fhd = str2func('cec17_func'); % Convert string to a function handle pointing to the CEC2017 benchmark entry function
func_nums = 30;          % Total number of benchmark functions available (informational here)
runs = 1;                % Number of independent runs per function

for i = 1:runs
    func_num = i;
    if func_num == 2     % Skip function #2 if encountered
        continue
    end
    for j = 1:runs
        [bestX, TotalBest, convergence_curve] = FSO(fhd, dim, lb, ub, maxIter, popSize, func_num);
    end
end

