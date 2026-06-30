%% Plots for 1 subject
total_trials = 9; % Number of recordings per movement
data = cell(1,total_trials);

% Load the subject data, found at ""DataCollection/trial[subject number]/p[movement number]%d.csv"
% [subject number] : goes from 1 to 4, each number is one different subject
% [movement number] : has 6 different possible values:
%   "" (nothing): Movement from origin to target A, without obstacle
%   1 : Movement from origin to target B, without obstacle
%   2 : Movement from origin to target C, without obstacle
%   3 : Movement from origin to target A, with obstacle
%   4 : Movement from origin to target B, with obstacle
%   5 : Movement from origin to target C, with obstacle
for i = 1:total_trials
    data{i} = load(sprintf("DataCollection/trial2/p%d.csv", i));
end

figure()
subplot(3,1,1)
% Process and plot the trials
firstPoint =zeros(1,total_trials); % Record first and last point, to obtain the variance across recordings
lastPoint =zeros(1,total_trials);
hold on
for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,1))
    title('X')
    
    firstPoint(i) = trialPosition(1,1);
    lastPoint(i) = trialPosition(end,1);
    %legend()
end
hold off
xFirstError = var(firstPoint);
xLastError = var(lastPoint);

subplot(3,1,2)
% Process and plot the trials
hold on
firstPoint =zeros(1,total_trials); % Record first and last point, to obtain the variance across recordings
lastPoint =zeros(1,total_trials);
for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,2))
    title('Y')
    firstPoint(i) = trialPosition(1,2);
    lastPoint(i) = trialPosition(end,2);
end
hold off
yFirstError = var(firstPoint);
yLastError = var(lastPoint);

subplot(3,1,3)
% Process and plot the trials
hold on
firstPoint =zeros(1,total_trials); % Record first and last point, to obtain the variance across recordings
lastPoint =zeros(1,total_trials);
for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,3))
    title('Z')
    firstPoint(i) = trialPosition(1,3);
    lastPoint(i) = trialPosition(end,3);
end
hold off
zFirstError = var(firstPoint);
zLastError = var(lastPoint);

fprintf('%%%%%%%%%%%%%%%%%% \n')
fprintf('X errors: %f %f \n', xFirstError, xLastError);
fprintf('Y errors: %f %f \n', yFirstError, yLastError);
fprintf('Z errors: %f %f \n', zFirstError, zLastError);



%% Plots for all subjects at the same time
total_trials = 9; % Number of recordings per movement

% Load the subject data, found at ""DataCollection/trial[subject number]/p[movement number]%d.csv"
% [subject number] : goes from 1 to 4, each number is one different subject
% [movement number] : has 6 different possible values:
%   "" (nothing): Movement from origin to target A, without obstacle
%   1 : Movement from origin to target B, without obstacle
%   2 : Movement from origin to target C, without obstacle
%   3 : Movement from origin to target A, with obstacle
%   4 : Movement from origin to target B, with obstacle
%   5 : Movement from origin to target C, with obstacle

data = cell(1,total_trials);
for i = 1:total_trials
    data{i} = load(sprintf("DataCollection/trial1/p%d.csv", i));
end

data2 = cell(1,total_trials);
for i = 1:total_trials
    data2{i} = load(sprintf("DataCollection/trial2/p%d.csv", i));
end

data3 = cell(1,total_trials);
for i = 1:total_trials
    data3{i} = load(sprintf("DataCollection/trial3/p%d.csv", i));
end

data4 = cell(1,total_trials);
for i = 1:total_trials
    data4{i} = load(sprintf("DataCollection/trial4/p%d.csv", i));
end

figure()
subplot(3,1,1)

hold on
for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,1),'r')
    title('X')

    trialTime = data2{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data2{i}(:,2:4); 
    plot(t_percent, trialPosition(:,1),'b')
    %legend()
    trialTime = data3{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data3{i}(:,2:4); 
    plot(t_percent, trialPosition(:,1),'g')

    trialTime = data4{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data4{i}(:,2:4); 
    plot(t_percent, trialPosition(:,1),'k')
end
hold off

subplot(3,1,2)

hold on

for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,2),'r')

    trialTime = data2{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data2{i}(:,2:4); 
    plot(t_percent, trialPosition(:,2),'b')
    %legend()
    trialTime = data3{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data3{i}(:,2:4);
    plot(t_percent, trialPosition(:,2),'g')

    trialTime = data4{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data4{i}(:,2:4); 
    plot(t_percent, trialPosition(:,2),'k')
    title('Y')

end
hold off


subplot(3,1,3)

hold on

for i = 1:total_trials
    trialTime = data{i}(:,1); % The first column is time
    t_percent = linspace(0, 100, length(trialTime)); % Normalize the time
    trialPosition = data{i}(:,2:4); % Columns 2 to 4 are X, Y and Z positions
    plot(t_percent, trialPosition(:,3),'r')

    trialTime = data2{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data2{i}(:,2:4); 
    plot(t_percent, trialPosition(:,3),'b')
    %legend()
    trialTime = data3{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data3{i}(:,2:4); 
    plot(t_percent, trialPosition(:,3),'g')

    trialTime = data4{i}(:,1);
    t_percent = linspace(0, 100, length(trialTime));
    trialPosition = data4{i}(:,2:4); 
    plot(t_percent, trialPosition(:,3),'k')
    title('Z')
end
hold off

