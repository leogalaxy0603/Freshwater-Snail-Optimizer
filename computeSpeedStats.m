% --------- Compute the average and standard deviation of speeds ---------
function [V_avg, sigma] = computeSpeedStats(velocities)
    % Compute the average speed
    V_avg = mean(vecnorm(velocities, 2, 2));  % Compute the Euclidean norm of each velocity vector, then take the mean
    % Compute the speed standard deviation
    sigma = std(vecnorm(velocities, 2, 2));   % Standard deviation of speeds
end
