% Verification script for the 2D heat transport model

% Initialize test parameters
test_grids = [50, 100, 200]; % Grid resolutions for testing
tolerance = 1e-6; % Tolerance for numerical error norms
validation_metrics = struct('MAE', [], 'RMSE', [], 'R2', []);

% Loop through test grids
for i = 1:length(test_grids)
    Nx = test_grids(i); % Grid resolution
    try
        % Run the main model with current resolution
        W = 16e3; % Domain width
        h = W / Nx; % Grid spacing
        run('./helmsdale_take_2.m'); % Main script
        
        % Check numerical error norms
        if (Errx > tolerance || Errz > tolerance)
            warning('Numerical error exceeds tolerance at Nx = %d', Nx);
        end

        % Validate against observed data
        temp_sim = interp1(zc, T, depth_obs);
        MAE = mean(abs(temp_obs - temp_sim));
        RMSE = sqrt(mean((temp_obs - temp_sim).^2));
        R2 = 1 - sum((temp_obs - temp_sim).^2) / sum((temp_obs - mean(temp_obs)).^2);

        % Store validation metrics
        validation_metrics(i).MAE = MAE;
        validation_metrics(i).RMSE = RMSE;
        validation_metrics(i).R2 = R2;

    catch ME
        fprintf('Error at Nx = %d: %s\n', Nx, ME.message);
    end
end

% Display validation summary
disp('Validation Metrics Summary:');
disp(validation_metrics);

% Plot validation metrics
figure;
plot(test_grids, [validation_metrics.MAE], 'o-', 'DisplayName', 'MAE');
hold on;
plot(test_grids, [validation_metrics.RMSE], 'x-', 'DisplayName', 'RMSE');
xlabel('Grid Resolution');
ylabel('Validation Metric');
legend;
title('Validation Metrics vs. Grid Resolution');
grid on;
