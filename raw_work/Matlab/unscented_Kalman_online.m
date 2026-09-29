%% Two-Stage UKF with Online Savitzky–Golay Velocity and Styled Plots
clear; close all; clc;

%% Load data
trialIdx = 5;
data = readmatrix('0 to 3_all_trials.xlsx');
cols = (trialIdx-1)*3 + (1:3);
time = data(:,cols(1));
pos  = data(:,cols(2));
valid = ~isnan(time) & ~isnan(pos);
time = time(valid); 
pos = pos(valid);
N = numel(time);
dt = median(diff(time));
if dt > 0.1
    time = time/1000; dt = dt/1000;
end
t_norm = linspace(0, 1, N);

%% UKF setup
nx = 3; %state dimensions
% parameters for dispersion
alpha = 1e-3; %spread of the points per state update
beta = 2; % gaussian distribution
kappa = 0; %secondary scaling
lambda = alpha^2*(nx + kappa) - nx;
wm = [lambda/(nx+lambda), repmat(1/(2*(nx+lambda)),1,2*nx)]; % weights for mean
wc = wm; % weights for covariance
wc(1) = wc(1) + (1 - alpha^2 + beta);
Q = diag([1e-5, 1e-5, 1e-4]); % process noise
R = diag([1e-4, 1e-4]);       % measurement noise

% Initialize state: [pos; vel; tau]
x = [pos(1); 0; pos(end)];
P = diag([1e-3, 1e-3, 1e-1]);

%% Parameters
H_steps = round(1000 / 1000 / dt);
future100_steps = round(0.1 / dt);
tau_est = nan(N,1);
pred_error = nan(N,1);
pred_horizon_ms = nan(N,1);

%% Video setup
fig = figure('Name', 'UKF Online Prediction');
vid = VideoWriter('UKF_online_prediction_visualized.mp4', 'MPEG-4');
vid.FrameRate = round(1/(dt*20));
open(vid);

%% Main filtering loop
for k = 2:N-1
    % === Online Savitzky–Golay Velocity ===
    frame_len = 9; poly_order = 3;
    M = min(frame_len, k);
    deg_eff = min(poly_order, M-1);
    t_win = time(k-M+1:k);
    p_win = pos(k-M+1:k);
    coef = polyfit(t_win - t_win(1), p_win, deg_eff);
    v_est = polyval(polyder(coef), t_win(end) - t_win(1));

    % Sigma points
    S = chol((nx+lambda)*P, 'lower');
    X = [x, x + S, x - S]; %three points, separed by S from the mean

    % Prediction
    X_pred = zeros(nx, 2*nx+1);
    for i = 1:size(X,2)
        xi = X(:,i);
        pos_p = xi(1) + xi(2)*dt; % next predicted pos is pos+vel*dt
        vel_p = xi(2) + 0.6*(xi(3) - xi(1))*dt; %the velocity will be higher the faster we are from the goal.
        tau_p = xi(3); % the goal stays the same
        X_pred(:,i) = [pos_p; vel_p; tau_p];
    end
    x_pred = X_pred * wm'; % the resulting prediction is the combination of all the samples
    P_pred = Q;
    for i = 1:size(X,2)
        dx = X_pred(:,i) - x_pred;
        P_pred = P_pred + wc(i) * (dx * dx'); % sample covariance
    end

    % Measurement update
    Z = X_pred(1:2,:);
    z_pred = Z * wm';
    P_zz = R;
    P_xz = zeros(nx,2);
    for i = 1:size(Z,2)
        dz = Z(:,i) - z_pred;
        dx = X_pred(:,i) - x_pred;
        P_zz = P_zz + wc(i) * (dz * dz');
        P_xz = P_xz + wc(i) * (dx * dz');
    end
    K = P_xz / P_zz; %update kalman gain
    z_meas = [pos(k); v_est];
    x = x_pred + K * (z_meas - z_pred);
    P = P_pred - K * P_zz * K';

    % Save goal estimate
    tau_k = x(3);
    tau_est(k) = tau_k;

    % Forecast future trajectory
    x_fore = x;
    pred_traj = nan(H_steps,1);
    for j = 1:H_steps
        x_fore(1) = x_fore(1) + x_fore(2)*dt;
        x_fore(2) = x_fore(2) + 0.6*(x_fore(3) - x_fore(1))*dt;
        pred_traj(j) = x_fore(1);
    end

    % Compute error and horizon
    future_len = min([future100_steps, H_steps, N - k]);
    pred_error(k) = sum(abs(pred_traj(1:future_len) - pos(k+1:k+future_len)));
    horizon_len = min([H_steps, N - k]);
    idxOK = find(abs(pred_traj(1:horizon_len) - pos(k+1:k+horizon_len)) <= deg2rad(1));
    pred_horizon_ms(k) = isempty(idxOK)*0 + (~isempty(idxOK))*idxOK(end)*dt*1000;

    % === Visualization ===
    clf;
    t_pred = t_norm(k)+(1:H_steps)*dt/time(end);
    t_pred = t_pred(1:min(H_steps,N-k));
    t_truth = t_norm(k+1:min(k+H_steps,N));

    plot(t_norm, pos, 'k--', 'LineWidth', 1.2); hold on;
    plot(t_norm(1:k), pos(1:k), 'b-', 'LineWidth', 2);
    plot(t_pred, pred_traj(1:numel(t_pred)), 'r-', 'LineWidth', 2);
    plot(t_truth, pos(k+1:min(k+H_steps,N)), 'g--', 'LineWidth', 2);
    plot(t_norm(k), pos(k), 'p', 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    legend('Full Ground Truth','Past Observed','Predicted','Actual Future','Current Point');
    xlabel('Normalized Time'); ylabel('Position (rad)');
    title(sprintf('Online Prediction: t = %.2fs', time(k)));
    ylim([min(pos)-0.1, max(pos)+0.1]);
    grid on;

    drawnow;
    writeVideo(vid,getframe(fig));
end
close(vid);

%% Post-animation plots

% Prediction error in 100 ms
figure;
plot(t_norm, pred_error, 'm-', 'LineWidth', 1.8);
xlabel('Normalized Time'); ylabel('Sum Abs Error (rad)');
title('Prediction Error in 100ms Horizon');
grid on;

% Prediction horizon
figure;
plot(t_norm, pred_horizon_ms, 'b-', 'LineWidth', 1.8);
xlabel('Normalized Time'); ylabel('Prediction Horizon (ms)');
title('Prediction Horizon with <1° Error');
grid on;

% Estimated τ vs true τ
figure;
plot(t_norm, tau_est, 'r', 'LineWidth', 2); hold on;
yline(pos(end), 'k--', 'LineWidth', 1.5);
xlabel('Normalized Time'); ylabel('Estimated Goal τ');
title('Online Estimated Goal vs True Goal');
legend('Estimated Goal','True Goal');
grid on;
