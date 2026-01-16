% --------- Collision checking function ---------
function [collisionCount, avgDistance] = checkCollisions(pop, gBest)
    % Set a threshold to determine whether particles are at the same position
    threshold = 0.5;  % A small threshold to judge “same position”
    [popSize, dim] = size(pop);
    
    % Compute Euclidean distances between particles
    % Use vectorized operations to compute distances among all particle pairs
    distanceMatrix = pdist2(pop, pop);  % Distance matrix; pdist2 computes pairwise distances
    collisionMatrix = distanceMatrix < threshold;  % Determine which particles are colliding
    
    % Zero out the diagonal to avoid comparing a particle with itself
    collisionMatrix = collisionMatrix - diag(diag(collisionMatrix));
    
    % Count the number of colliding pairs
    collisionCount = sum(collisionMatrix(:)) / 2;  % Symmetric matrix; divide by 2 to avoid double-counting
    
    % Compute the average distance from colliding particles to the global best gBest
    collisionParticles = find(any(collisionMatrix, 2));  % Indices of particles involved in collisions
    
    if isempty(collisionParticles)
        avgDistance = 0;  % If there are no colliding particles, return 0
    else
        % Compute distances from colliding particles to the global best
        distancesToGBest = vecnorm(pop(collisionParticles, :) - gBest, 2, 2);
        avgDistance = mean(distancesToGBest);  % Average distance
    end
end
