function [bestSolution, bestFitness, convergence_curve] = FSO(fhd, dim, lb, ub, maxIter, popSize,varargin)
% Freshwater Snail Optimization Algorithm (FSO)
% This function implements single-objective continuous optimization inspired by
% freshwater snails moving in water flow.
% Usage:
%    [bestSolution, bestFitness] = snail_optimize(func, dim, lb, ub, maxIter, popSize)
% Parameters:
%    func     - Objective function handle (@func) or function name string.
%               The objective should accept a 1×dim vector and return a scalar fitness (smaller is better).
%    dim      - Problem dimension (number of decision variables).
%    lb       - Lower bounds of decision variables (scalar or 1×dim vector).
%    ub       - Upper bounds of decision variables (scalar or 1×dim vector).
%    maxIter  - Maximum number of iterations.
%    popSize  - Population size (number of snails).
% Outputs:
%    bestSolution - Best solution found (1×dim vector).
%    bestFitness  - Objective value of the best solution.
%
% Expand bounds to row vectors of length dim
lb = lb(:)';  ub = ub(:)';
if numel(lb) == 1
    lb = lb * ones(1, dim);
end
if numel(ub) == 1
    ub = ub * ones(1, dim);
end
if length(lb) ~= dim || length(ub) ~= dim
    error('The lower bound lb and upper bound ub should be scalars or vectors of length dim.');
end
beta = 7;          % Exponential decay rate controlling how A decreases over iterations
% Initialize snail population with random positions (uniform over the search space)
snailPos =  initialization(popSize, dim, ub, lb);
% Compute initial fitness of the population
snailFitness = zeros(popSize, 1);
for i = 1:popSize
    snailFitness(i) = feval(fhd,snailPos(i, :)',varargin{:});
end
% Initialize each individual's personal best with its initial state
snailPbestPos = snailPos;
snailPbestFit = snailFitness;
% Find the initial global best
[bestFitness, bestIndex] = min(snailFitness);
bestSolution = snailPos(bestIndex, :);

initialExploreRate = 0.7;
finalExploreRate   = 0.1;
P_collisionCount = 0;
P_avgDistance =  0;

% Main optimization loop
for iter = 1:maxIter

    [V_avg, sigma] = computeSpeedStats(bestSolution);  % Compute average speed and standard deviation for the current generation
    gamma = (4 / pi) * (1 + 0.4 * (sigma / V_avg)^2); % Compute alpha
    collisionCount_gen(iter,:) = checkCollisions(snailPos, bestSolution);
    [collisionCount, avgDistance] = checkCollisions(snailPos, bestSolution);

    if collisionCount < popSize
        lambda = 1/(2*gamma*avgDistance*collisionCount);
        lambda_gen(iter,:) = lambda;

        P_collisionCount = collisionCount;
        P_avgDistance = avgDistance;
        %
    elseif P_collisionCount == 0 && P_avgDistance == 0
        lambda = 0;
        lambda_gen(iter,:) = lambda;
    else
        % When collisions almost occur, use this mode for a free step-length strategy
        collisionCount = popSize - P_collisionCount;
        lambda = mean(bestSolution)/(2*P_avgDistance*V_avg*collisionCount);
        lambda_gen(iter,:) = lambda;

    end
    snailPos = snailPos+snailPos.*lambda;

    % Dynamically compute the number of explorers in this generation
    exploreRate = initialExploreRate - (initialExploreRate - finalExploreRate) * (iter / maxIter);
    explorerCount = round(exploreRate * popSize);
    if explorerCount < 1
        explorerCount = 1;
    elseif explorerCount > popSize
        explorerCount = popSize;
    end

    % Sort by current fitness (high to low) to determine explorer indices
    [~, sortIndex] = sort(snailFitness, 'descend'); 
    explorerIdx = sortIndex(1:explorerCount);
    attachedIdx = sortIndex(explorerCount+1:end);

    % Save current global best position for guidance (previous generation's global best)
    gbestPos = bestSolution;
    % Preallocate new position matrix
    newPos = zeros(popSize, dim);

    % --- Explorers: drifting (global exploration) ---
    for idx = explorerIdx'
        A0 = rand();             % Initial amplitude of the influence factor (upper bound of the step-size coefficient)
        A = A0 * exp(-beta * iter / maxIter);        % Exponentially decaying step-size coefficient A (simulate decreasing buoyancy)
        A_gen(iter)=A;
        % Randomly reset to a new position within the search space
        if rand<0.5
            sine_wave_influence = sin(2 * pi * norm(gbestPos - snailPbestPos(idx, :)) / (ub(:,1) - lb(:,1)));  % Fluctuation based on positions
        else
            sine_wave_influence = cos(2 * pi * norm(gbestPos - snailPbestPos(idx, :)) / (ub(:,1) - ub(:,1)));  % Fluctuation based on positions
        end
        newPos(idx, :) = A*(gbestPos-snailPos(idx, :))+ sine_wave_influence .*(snailPbestPos(idx, :)-snailPos(idx,:));
    end

    % --- Exploiters: adhesion crawling (local exploitation) ---
    % Compute the local step-size factor and global guidance factor for this iteration
    alpha = 1 - iter / maxIter;         % Coefficient decreasing with iterations (from 1 down to 0)

    localFactor  = 0.1 * alpha;         % Local random step-size coefficient (starts at 0.1, decreases to 0)
    globalFactor = 0.2 * (1 - alpha);   % Global guidance coefficient (starts at 0, increases to 0.2)
    for idx = attachedIdx'
        % Small random perturbation around the individual's personal best
        % Also add a small offset toward the current global best
        randStep = (rand(1, dim) * 2 - 1) .* (ub - lb); 
        newPos(idx, :) = snailPbestPos(idx, :) + localFactor .* randStep + globalFactor .* (gbestPos - snailPbestPos(idx, :));
        % Boundary control: ensure the new position stays within [lb, ub]
        newPos(idx, :) = max(newPos(idx, :), lb);
        newPos(idx, :) = min(newPos(idx, :), ub);
    end

    % Apply boundary control to explorers as well (keep all dimensions within valid ranges)
    for idx = explorerIdx'
        newPos(idx, :) = max(newPos(idx, :), lb);
        newPos(idx, :) = min(newPos(idx, :), ub);
    end

    % Evaluate the new generation and update personal/global bests
    for i = 1:popSize
        % Compute the new fitness of the current individual
        fval = feval(fhd,newPos(i, :)',varargin{:});
        snailFitness(i) = fval;
        % Update personal best
        if fval < snailPbestFit(i)
            snailPbestFit(i) = fval;
            snailPbestPos(i, :) = newPos(i, :);
        end
        % Update global best
        if fval < bestFitness
            bestFitness = fval;
            bestSolution = newPos(i, :);
        end
    end
    if mod(iter, 20) == 0
        disp(['Gen: ' num2str(iter) ', BestFitness: ' num2str(bestFitness)]);
    end
    % Update population positions to the new positions
    snailPos = newPos;
    convergence_curve(iter)=bestFitness;
end

end
