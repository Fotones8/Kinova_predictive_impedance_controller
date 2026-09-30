% data {mode, movement}
modes = 6;
movements = 5;
trials = 5;
subjects = 15;


data = cell(modes,trials,movements,subjects);
for j = 1:modes
    for i = 1:movements
        for z = 1:trials
            for s = 1:subjects
            data{j,i,z,s} = readmatrix(sprintf("subject%d/subject%dmode%dmov%d%d.csv", s+30, s+30, j,z,i));
        %fprintf('%d %d\n', j, i)
            end
        end
    end
end

data_unclipped = data;


%% Clip data based on total velocity threshold
data = data_unclipped;

data_clipped = cell(modes,trials,movements,subjects);
for j = 1:modes
    for i = 1:trials
        for z = 1:movements
            for s = 1:subjects
            current_data = data{j,i,z,s};
            fprintf('%d %d %d %d\n', s, j, i,z)
            
            % Compute total velocity magnitude
            total_velocity = sqrt(current_data(:,8).^2 + current_data(:,9).^2 + current_data(:,10).^2);
            
            % Threshold = 10% of max velocity
            threshold = 0.10 * max(total_velocity);
            
            % Find the peak velocity index
            [~, peak_idx] = max(total_velocity);
            
            % Walk backward from peak to find start_idx
            start_idx = peak_idx;
            while start_idx > 1 && total_velocity(start_idx - 1) >= threshold
                start_idx = start_idx - 1;
            end
            
            % Walk forward from peak to find end_idx
            end_idx = peak_idx;
            n = length(total_velocity);
            while end_idx < n && total_velocity(end_idx + 1) >= threshold
                end_idx = end_idx + 1;
            end
            
            data_clipped{j,i,z,s} = current_data(start_idx:end_idx, :);
            end
        end
    end
end


j=1; i=1; z=1;
figure
subplot(2,1,1)
plot(sqrt(data{2,2,2,1}(:,8).^2+data{2,2,2,1}(:,9).^2+data{2,2,2,1}(:,10).^2))
title('Original')
subplot(2,1,2)
plot(sqrt(data_clipped{2,2,2,1}(:,8).^2+data_clipped{2,2,2,1}(:,9).^2+data_clipped{2,2,2,1}(:,10).^2))
title('Clipped')
data = data_clipped;


%%
figure()
hold on

r = 0.02;   % radius (diameter 0.04)
[Xs, Ys, Zs] = sphere(20);
center = [0.51, 0.16, 0.43];

surf(r*Xs + center(1), r*Ys + center(2), r*Zs + center(3), ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none')

r = 0.02;
[Xs, Ys, Zs] = sphere(20);
desaturationFactor = 0.35;   % 0 = fully gray, 1 = fully original color

for i = 1:movements
    traj = data{1,i,2,1};
    h = plot3(traj(:,2), traj(:,3), traj(:,4));
    lineColor = h.Color;

    disp(i)
    disp(traj(1,2))
    disp(traj(1,3))
    disp(traj(1,4))
    
    hsvColor = rgb2hsv(lineColor);
    hsvColor(2) = hsvColor(2) * desaturationFactor;   % reduce saturation
    ballColor = hsv2rgb(hsvColor);
    
    endPoint = traj(end, 2:4);
    surf(r*Xs + endPoint(1), r*Ys + endPoint(2), r*Zs + endPoint(3), ...
        'FaceColor', ballColor, 'EdgeColor', 'none')
end

hold off
axis equal
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')

%% Human effort
%groups = [1,2,3,4,5];
grasp_force_mean =[];
subj = [];
controller  = [];
goal = [];
for j = 1:modes
    for i = 1:movements
        data_all_trials = [];
        for t = 1:trials
            for s = 1:subjects
        %data_all_trials = [data{j,i,1}(1:end-1,36);data{j,i,2}(1:end-1,36);data{j,i,3}(1:end-1,36);data{j,i,4}(1:end-1,36);data{j,i,5}(1:end-1,36)];
                grasp_force_mean = [grasp_force_mean; mean(data{j,i,t,s}(1:end-1,36))];
                subj = [subj; s];
                controller = [controller; j];
                goal = [goal; i];
         
            end
        end
    end
end


[p, tbl, stats] = anovan(grasp_force_mean, {subj, controller, goal}, ...
    'model', 'interaction', ...
    'varnames', {'Subject','Controller','Target'})

% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
[results_controller, m, h, gnames] = multcompare(stats, 'Dimension', 2)

groupMeans = m(:,1);
groupSE = m(:,2);

% 95% CI half-width using the t-distribution (matches the CI convention multcompare uses)
alpha = 0.05;
critVal = tinv(1 - alpha/2, stats.dfe);
ciHalfWidth = critVal * groupSE;

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};

% Reorder: original index order 1,6,2,3,4,5
newOrder = [1, 6, 2, 3, 4, 5];
groupMeans = groupMeans(newOrder);
ciHalfWidth = ciHalfWidth(newOrder);
customLabels = customLabels(newOrder);


figure
bar(groupMeans)
hold on
errorbar(1:length(groupMeans), groupMeans, ciHalfWidth, 'k', 'LineStyle', 'none', 'LineWidth', 1.2)
hold off

set(gca, 'XTick', 1:length(groupMeans), 'XTickLabel', customLabels,'FontSize', 14)
ylabel('Mean grasp force (N)','FontSize', 16)
%title('Mean grasp force by Controller (95% CI)')


% Pairwise comparisons for Goal (dimension 3)
figure
results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);


customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];

controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

figure
hold on
boxchart(controllerCat, grasp_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

jitterWidth = 0.3;
xNumeric = double(controllerCat) + (rand(length(controllerCat),1)-0.5)*jitterWidth;
scatter(xNumeric, grasp_force_mean, 4, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.3)

hold off
set(gca, 'FontSize', 14)
ylabel('Mean grasp force (N)', 'FontSize', 16)
ylim([-2,40])

%% New graph human effort
customLabels = {'FA-1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'FA-2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

% Collapse to one mean per subject per controller (averaging over movements/trials)
subjControllerPairs = unique([subj, controller], 'rows');
n_pairs = size(subjControllerPairs, 1);

subjectMeans = zeros(n_pairs, 1);
pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));

for k = 1:n_pairs
    thisSubj = subjControllerPairs(k,1);
    thisController = subjControllerPairs(k,2);
    mask = (subj == thisSubj) & (controller == thisController);
    subjectMeans(k) = mean(grasp_force_mean(mask));
end

% Box plot (unchanged) + one point per subject
figure
hold on
boxchart(controllerCat, grasp_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5],'LineWidth', 1.5)

jitterWidth = 0.3;
xNumeric = double(pairControllerCat) + (rand(n_pairs,1)-0.5)*jitterWidth;
scatter(xNumeric, subjectMeans, 24, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6)

% Compact letter display: groups sharing a letter are NOT significantly different
letterLabels = {'A', 'B', 'C', 'CD', 'CD', 'D'};   % matches display order: FA1, FA2, ProMP var, UKF var, ProMP fixed, UKF fixed

whiskerTop = zeros(1,6);
for g = 1:6
    groupVals = grasp_force_mean(double(controllerCat) == g);
    q1 = prctile(groupVals, 25);
    q3 = prctile(groupVals, 75);
    iqrVal = q3 - q1;
    upperFence = q3 + 1.5*iqrVal;
    whiskerTop(g) = max(groupVals(groupVals <= upperFence));
end


yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(grasp_force_mean(double(controllerCat) == g));
end
%letterY = yMaxPerGroup + 2;   % offset above each box's max point
letterY = whiskerTop + 2;
letterY(2) = 25;

for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold')
end

% % Explicit bracket for the one notable within-group difference (ProMP var vs UKF fixed, p=0.0147)
% x1 = 3;  % ProMP var
% x2 = 6;  % UKF fixed
% %yBracket = max(yMaxPerGroup) + 6;
% yBracket = 27;
% 
% plot([x1 x1 x2 x2], [yBracket-0.5, yBracket, yBracket, yBracket-0.5], 'k-', 'LineWidth', 1)
% text(mean([x1 x2]), yBracket + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold')

hold off
set(gca, 'FontSize', 16)
ylabel('Grasping force (N)', 'FontSize', 18)
ylim([-5, 30])   % extended to fit the letters and bracket

%% %%%%%%%%%%% with std dev

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

% Collapse to one mean (+ std) per subject per controller
subjControllerPairs = unique([subj, controller], 'rows');
n_pairs = size(subjControllerPairs, 1);

subjectMeans = zeros(n_pairs, 1);
subjectStds = zeros(n_pairs, 1);
pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));

for k = 1:n_pairs
    thisSubj = subjControllerPairs(k,1);
    thisController = subjControllerPairs(k,2);
    mask = (subj == thisSubj) & (controller == thisController);
    subjectMeans(k) = mean(grasp_force_mean(mask));
    subjectStds(k) = std(grasp_force_mean(mask));
end

% Box plot + one point per subject with std error bars
figure
hold on
boxchart(controllerCat, grasp_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

jitterWidth = 0.3;
xNumeric = double(pairControllerCat) + (rand(n_pairs,1)-0.5)*jitterWidth;

errorbar(xNumeric, subjectMeans, subjectStds, 'LineStyle', 'none', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, 'CapSize', 0)
scatter(xNumeric, subjectMeans, 12, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6)

% Compact letter display: groups sharing a letter are NOT significantly different
letterLabels = {'A', 'B', 'C', 'CD', 'CD', 'D'};   % matches display order: FA1, FA2, ProMP var, UKF var, ProMP fixed, UKF fixed

yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(grasp_force_mean(double(controllerCat) == g));
end
letterY = yMaxPerGroup + 2;   % offset above each box's max point

for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold')
end

% Explicit bracket for the one notable within-group difference (ProMP var vs UKF fixed, p=0.0147)
x1 = 3;  % ProMP var
x2 = 6;  % UKF fixed
yBracket = max(yMaxPerGroup) + 6;

plot([x1 x1 x2 x2], [yBracket-0.5, yBracket, yBracket, yBracket-0.5], 'k-', 'LineWidth', 1)
text(mean([x1 x2]), yBracket + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold')

hold off
set(gca, 'FontSize', 14)
ylabel('Mean grasp force (N)', 'FontSize', 16)
ylim([-5, 60])   % extended to fit the letters and bracket

%%
customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

% Collapse to one mean (+ std) per subject per controller
subjControllerPairs = unique([subj, controller], 'rows');
n_pairs = size(subjControllerPairs, 1);
subjectMeans = zeros(n_pairs, 1);
subjectStds = zeros(n_pairs, 1);
pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));

for k = 1:n_pairs
    thisSubj = subjControllerPairs(k,1);
    thisController = subjControllerPairs(k,2);
    mask = (subj == thisSubj) & (controller == thisController);
    subjectMeans(k) = mean(grasp_force_mean(mask));
    subjectStds(k) = std(grasp_force_mean(mask));
end

% Regular, evenly-spaced x-offsets per subject within each controller group
spreadWidth = 0.45;   % total horizontal spread of the 15 points within a group's box
xNumeric = zeros(n_pairs, 1);

controllerCodes = double(pairControllerCat);
uniqueCodes = unique(controllerCodes);

for g = uniqueCodes'
    groupIdx = find(controllerCodes == g);
    [~, sortOrder] = sort(subjControllerPairs(groupIdx, 1));   % order by subject ID
    orderedIdx = groupIdx(sortOrder);
    
    nSubjInGroup = length(orderedIdx);
    offsets = linspace(-spreadWidth/2, spreadWidth/2, nSubjInGroup);
    
    xNumeric(orderedIdx) = g + offsets;
end

% Box plot + one point per subject with std error bars
figure
hold on
boxchart(controllerCat, grasp_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

errorbar(xNumeric, subjectMeans, subjectStds, 'LineStyle', 'none', 'Color', [0 0.2 0.55], 'LineWidth', 0.8, 'CapSize', 3)
scatter(xNumeric, subjectMeans, 12, [0 0.2 0.55], 'filled', 'MarkerFaceAlpha', 0.6)

% Compact letter display: groups sharing a letter are NOT significantly different
letterLabels = {'A', 'B', 'C', 'CD', 'CD', 'D'};   % matches display order: FA1, FA2, ProMP var, UKF var, ProMP fixed, UKF fixed
yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(grasp_force_mean(double(controllerCat) == g));
end
letterY = yMaxPerGroup + 2;
for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold')
end

% Explicit bracket for ProMP var vs UKF fixed (p=0.0147)
x1 = 3; x2 = 6;
yBracket = max(yMaxPerGroup) + 6;
plot([x1 x1 x2 x2], [yBracket-0.5, yBracket, yBracket, yBracket-0.5], 'k-', 'LineWidth', 1)
text(mean([x1 x2]), yBracket + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold')

hold off
set(gca, 'FontSize', 14)
ylabel('Mean grasp force (N)', 'FontSize', 16)
ylim([-5, 60])

%% Human effort graph

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

figure
hold on
boxchart(controllerCat, grasp_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

jitterWidth = 0.3;
xNumeric = double(controllerCat) + (rand(length(controllerCat),1)-0.5)*jitterWidth;
scatter(xNumeric, grasp_force_mean, 4, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.3)

% Compact letter display: groups sharing a letter are NOT significantly different
letterLabels = {'A', 'B', 'C', 'CD', 'CD', 'D'};   % matches display order: FA1, FA2, ProMP var, UKF var, ProMP fixed, UKF fixed

yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(grasp_force_mean(double(controllerCat) == g));
end
letterY = yMaxPerGroup + 2;   % offset above each box's max point

for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold')
end

% Explicit bracket for the one notable within-group difference (ProMP var vs UKF fixed, p=0.0147)
x1 = 3;  % ProMP var
x2 = 6;  % UKF fixed
yBracket = max(yMaxPerGroup) + 6;

plot([x1 x1 x2 x2], [yBracket-0.5, yBracket, yBracket, yBracket-0.5], 'k-', 'LineWidth', 1)
text(mean([x1 x2]), yBracket + 0.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold')

hold off
set(gca, 'FontSize', 14)
ylabel('Mean grasp force (N)', 'FontSize', 16)
ylim([-5, 60])   % extended to fit the letters and bracket

%%
% function results = checkNormalityByGroup(data, groupVar, groupName)
%     levels = unique(groupVar);
%     results = table();
%     for c = 1:length(levels)
%         groupData = data(groupVar == levels(c));
%         [H, pValue, W] = swtest(groupData);
%         results.(groupName)(c) = levels(c);
%         results.N(c) = length(groupData);
%         results.W(c) = W;
%         results.pValue(c) = pValue;
%         results.RejectNormality(c) = H;
%     end
% end
% 
% resultsByController = checkNormalityByGroup(grasp_force_mean, controller, 'Controller');
% disp(resultsByController)
% resultsBySubject = checkNormalityByGroup(grasp_force_mean, subj, 'Subject');
% disp(resultsBySubject)
% resultsByGoal = checkNormalityByGroup(grasp_force_mean, goal, 'Target');
% disp(resultsByGoal)

%% Interaction force
%groups = [1,2,3,4,5,6];
%movements=5;

Interaction_force_mean = [];
%Interaction_force_max = zeros(modes,movements);

Interaction_force_meanX = [];
Interaction_force_meanY = [];
Interaction_force_meanZ = [];

subj = [];
controller = [];
goal = [];
baseline_Z = -27;

for j = 1:modes
    for i = 1:movements
        data_all_trials = [];
        for t = 1:trials
            for s = 1:subjects
                Interaction_force_mean = [Interaction_force_mean; mean(sqrt(data{j,i,t,s}(1:end-1,14).^2+data{j,i,t,s}(1:end-1,15).^2+(data{j,i,t,s}(1:end-1,16)-baseline_Z).^2))];        
                Interaction_force_meanX = [Interaction_force_meanX; mean(sqrt(data{j,i,t,s}(1:end-1,14).^2))];
                Interaction_force_meanY = [Interaction_force_meanY; mean(sqrt(data{j,i,t,s}(1:end-1,15).^2))];
                Interaction_force_meanZ = [Interaction_force_meanZ; mean(sqrt((data{j,i,t,s}(1:end-1,16)-baseline_Z).^2))];
                subj = [subj; s];
                controller = [controller; j];

                goal = [goal; i];
         
            end
        end
    end
end


[p, tbl, stats] = anovan(Interaction_force_mean, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
[results_controller, m, h, gnames] = multcompare(stats, 'Dimension', 2);


groupMeans = m(:,1);
groupSE = m(:,2);

% 95% CI half-width using the t-distribution (matches the CI convention multcompare uses)
alpha = 0.05;
critVal = tinv(1 - alpha/2, stats.dfe);
ciHalfWidth = critVal * groupSE;

customLabels = {'Force-amplification 1', 'ProMP var', 'Kalman var', 'ProMP fixed', 'Kalman fixed', 'Force-amplification 2'};

% Reorder: original index order 1,6,2,3,4,5
newOrder = [1, 6, 2, 3, 4, 5];
groupMeans = groupMeans(newOrder);
ciHalfWidth = ciHalfWidth(newOrder);
customLabels = customLabels(newOrder);

figure
bar(groupMeans)
hold on
errorbar(1:length(groupMeans), groupMeans, ciHalfWidth, 'k', 'LineStyle', 'none', 'LineWidth', 1.2)
hold off



set(gca, 'XTick', 1:length(groupMeans), 'XTickLabel', customLabels,'FontSize', 16)
ylabel('Mean interaction force (N)','FontSize', 16)
%title('Mean interaction force by Controller (95% CI)')




% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
%figure
%results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
%figure
%results_interaction = multcompare(stats, 'Dimension', [1 2]);

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

figure
hold on
boxchart(controllerCat, Interaction_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

jitterWidth = 0.3;
xNumeric = double(controllerCat) + (rand(length(controllerCat),1)-0.5)*jitterWidth;
scatter(xNumeric, Interaction_force_mean, 4, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.3)

% Compact letter display: FA1 & FA2 share a letter, middle four share a different letter
% Display order: FA1, FA2, ProMP var, Kalman var, ProMP fixed, Kalman fixed
letterLabels = {'A', 'A', 'B', 'B', 'B', 'B'};

yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(Interaction_force_mean(double(controllerCat) == g));
end
letterY = min(yMaxPerGroup + 0.05*range(Interaction_force_mean), 40);   % adjust offset to your data's scale

for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold')
end

hold off
set(gca, 'FontSize', 14)
ylabel('Mean interaction force (N)', 'FontSize', 16)
%ylim([min(Interaction_force_mean), max(Interaction_force_mean) + 0.15*range(Interaction_force_mean)])
ylim([0,45])
%% With means per subject INTERACTION FORCE

customLabels = {'FA-1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'FA-2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

% Collapse to one mean (+ std) per subject per controller
subjControllerPairs = unique([subj, controller], 'rows');
n_pairs = size(subjControllerPairs, 1);
subjectMeans = zeros(n_pairs, 1);
subjectStds = zeros(n_pairs, 1);
pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));

for k = 1:n_pairs
    thisSubj = subjControllerPairs(k,1);
    thisController = subjControllerPairs(k,2);
    mask = (subj == thisSubj) & (controller == thisController);
    subjectMeans(k) = mean(Interaction_force_mean(mask));
    subjectStds(k) = std(Interaction_force_mean(mask));
end

% Regular, evenly-spaced x-offsets per subject within each controller group
spreadWidth = 0.45;
xNumeric = zeros(n_pairs, 1);
controllerCodes = double(pairControllerCat);
uniqueCodes = unique(controllerCodes);

for g = uniqueCodes'
    groupIdx = find(controllerCodes == g);
    [~, sortOrder] = sort(subjControllerPairs(groupIdx, 1));
    orderedIdx = groupIdx(sortOrder);
    nSubjInGroup = length(orderedIdx);
    offsets = linspace(-spreadWidth/2, spreadWidth/2, nSubjInGroup);
    xNumeric(orderedIdx) = g + offsets;
end

% Box plot + one jittered point per subject
figure
hold on
boxchart(controllerCat, Interaction_force_mean, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5], 'LineWidth', 1.5)

jitterWidth = 0.3;
xNumeric = double(pairControllerCat) + (rand(n_pairs,1)-0.5)*jitterWidth;
scatter(xNumeric, subjectMeans, 24, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6)

% Compact letter display, placed just above each group's upper whisker
letterLabels = {'A', 'A', 'B', 'B', 'B', 'B'};

whiskerTop = zeros(1,6);
for g = 1:6
    groupVals = Interaction_force_mean(double(controllerCat) == g);
    q1 = prctile(groupVals, 25);
    q3 = prctile(groupVals, 75);
    iqrVal = q3 - q1;
    upperFence = q3 + 1.5*iqrVal;
    whiskerTop(g) = max(groupVals(groupVals <= upperFence));
end

letterY = min(whiskerTop + 1, 40);
for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold')
end

hold off
set(gca, 'FontSize', 16)
ylabel('Interaction force (N)', 'FontSize', 18)
ylim([0,25])


%%
[p, tbl, stats] = anovan(Interaction_force_meanX, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
[results_controller, m, h, gnames] = multcompare(stats, 'Dimension', 2);

figure()
groupMeans = m(:,1);
bar(groupMeans)

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
set(gca, 'XTick', 1:length(groupMeans), 'XTickLabel', customLabels,'FontSize', 14)
% Pairwise comparisons for Subject (dimension 1)
%figure
%results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
%figure
%results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
%figure
%results_interaction = multcompare(stats, 'Dimension', [1 2]);
%%
[p, tbl, stats] = anovan(Interaction_force_meanY, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
[results_controller, m, h, gnames] = multcompare(stats, 'Dimension', 2);

figure()
groupMeans = m(:,1);
bar(groupMeans)

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
set(gca, 'XTick', 1:length(groupMeans), 'XTickLabel', customLabels,'FontSize', 14)
% Pairwise comparisons for Subject (dimension 1)
%figure
%results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
%figure
%results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
%figure
%results_interaction = multcompare(stats, 'Dimension', [1 2]);
%%
[p, tbl, stats] = anovan(Interaction_force_meanZ, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
%figure
%results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
%figure
%results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
%figure
%results_interaction = multcompare(stats, 'Dimension', [1 2]);
%%
%% Build long-format data: one row per trial per axis
n = length(controller);

allValues = [Interaction_force_meanX; Interaction_force_meanY; Interaction_force_meanZ];
axisLabels = [repmat("X", n, 1); repmat("Y", n, 1); repmat("Z", n, 1)];
controllerAll = [controller; controller; controller];

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];

controllerCat = categorical(controllerAll, [1 6 2 3 4 5], customLabels(newOrder));
axisCat = categorical(axisLabels, {'X','Y','Z'});


% Box plot
figure
hold on
colororder([0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3])
boxchart(controllerCat, allValues, 'GroupByColor', axisCat, 'MarkerStyle', 'none')

axisNamesList = {'X','Y','Z'};
axisOffsets = containers.Map({'X','Y','Z'}, {-0.25, 0, 0.25});
axisColors = [0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3];
jitterWidth = 0.15;

for a = 1:3
    thisAxisName = axisNamesList{a};
    mask = (axisLabels == thisAxisName);
    
    xNumeric = double(controllerCat(mask)) + axisOffsets(thisAxisName) + (rand(sum(mask),1)-0.5)*jitterWidth;
    
    scatter(xNumeric, allValues(mask), 4, axisColors(a,:), 'filled', 'MarkerFaceAlpha', 0.3)
end

hold off
set(gca, 'FontSize', 14)
ylabel('Interaction force (N)', 'FontSize', 16)
ylim([0, 30])
legend({'X','Y','Z'}, 'Location', 'best')

%% INTERACTION FORCES 3 AXIS, COMPLETE SCATTER
n = length(controller);
allValues = [Interaction_force_meanX; Interaction_force_meanY; Interaction_force_meanZ];
axisLabels = [repmat("X", n, 1); repmat("Y", n, 1); repmat("Z", n, 1)];
controllerAll = [controller; controller; controller];

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controllerAll, [1 6 2 3 4 5], customLabels(newOrder));
axisCat = categorical(axisLabels, {'X','Y','Z'});

figure
hold on
colororder([0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3])
boxchart(controllerCat, allValues, 'GroupByColor', axisCat, 'MarkerStyle', 'none')

axisNamesList = {'X','Y','Z'};
axisOffsets = containers.Map({'X','Y','Z'}, {-0.25, 0, 0.25});
axisColors = [0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3];
jitterWidth = 0.15;

for a = 1:3
    thisAxisName = axisNamesList{a};
    mask = (axisLabels == thisAxisName);
    xNumeric = double(controllerCat(mask)) + axisOffsets(thisAxisName) + (rand(sum(mask),1)-0.5)*jitterWidth;
    scatter(xNumeric, allValues(mask), 4, axisColors(a,:), 'filled', 'MarkerFaceAlpha', 0.3)
end

% Compact letter display, per axis
% Display order: FA1, FA2, ProMP var, UKF var, ProMP fixed, UKF fixed
lettersX = {'A', 'A', 'B', 'C', 'B', 'C'};
lettersY = {'A', 'A', 'B', 'C', 'B', 'D'};
lettersZ = {'AB', 'A', 'AB', 'B', 'AB', 'C'};

axisLetterSets = {lettersX, lettersY, lettersZ};
axisData = {Interaction_force_meanX, Interaction_force_meanY, Interaction_force_meanZ};

for a = 1:3
    letters = axisLetterSets{a};
    valsThisAxis = axisData{a};
    xOffset = axisOffsets(axisNamesList{a});
    
    for g = 1:6
        groupVals = valsThisAxis(double(controllerCat(1:n)) == g);  % controller labels repeat every n rows
        yMaxGroup = max(groupVals);
        
        text(g + xOffset, min(yMaxGroup + 0.5, 32), letters{g}, ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', axisColors(a,:))
    end
end

hold off
set(gca, 'FontSize', 14)
ylabel('Interaction force (N)', 'FontSize', 16)
ylim([0, 33])
legend({'X','Y','Z'}, 'Location', 'best')

%% INTERACTION FORCES 3 AXIS SUBJECT MEAN
n = length(controller);
allValues = [Interaction_force_meanX; Interaction_force_meanY; Interaction_force_meanZ];
axisLabels = [repmat("X", n, 1); repmat("Y", n, 1); repmat("Z", n, 1)];
controllerAll = [controller; controller; controller];
customLabels = {'FA-1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'FA-2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controllerAll, [1 6 2 3 4 5], customLabels(newOrder));
axisCat = categorical(axisLabels, {'X','Y','Z'});

figure
hold on
colororder([0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3])
boxchart(controllerCat, allValues, 'GroupByColor', axisCat, 'MarkerStyle', 'none', 'LineWidth', 1.5)

axisNamesList = {'X','Y','Z'};
axisOffsets = containers.Map({'X','Y','Z'}, {-0.22, 0, 0.22});
axisColors = [0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3];
axisData = {Interaction_force_meanX, Interaction_force_meanY, Interaction_force_meanZ};
spreadWidth = 0.18;

for a = 1:3
    thisAxisName = axisNamesList{a};
    valsThisAxis = axisData{a};
    xOffset = axisOffsets(thisAxisName);
    subjControllerPairs = unique([subj, controller], 'rows');
    n_pairs = size(subjControllerPairs, 1);
    subjectMeans = zeros(n_pairs, 1);
    subjectStds = zeros(n_pairs, 1);
    pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));
    
    for k = 1:n_pairs
        thisSubj = subjControllerPairs(k,1);
        thisController = subjControllerPairs(k,2);
        mask = (subj == thisSubj) & (controller == thisController);
        subjectMeans(k) = mean(valsThisAxis(mask));
        subjectStds(k) = std(valsThisAxis(mask));
    end
    
    jitterWidth = 0.1;
    xNumeric = double(pairControllerCat) + axisOffsets(thisAxisName) + (rand(n_pairs,1)-0.5)*jitterWidth;
    scatter(xNumeric, subjectMeans, 24, axisColors(a,:), 'filled', 'MarkerFaceAlpha', 0.6)
end

% Compact letter display, per axis — placed just above each group's whisker
lettersX = {'A', 'A', 'B', 'C', 'B', 'C'};
lettersY = {'A', 'A', 'B', 'C', 'B', 'D'};
lettersZ = {'AB', 'A', 'AB', 'B', 'AB', 'C'};
axisLetterSets = {lettersX, lettersY, lettersZ};

for a = 1:3
    letters = axisLetterSets{a};
    valsThisAxis = axisData{a};
    xOffset = axisOffsets(axisNamesList{a});
    for g = 1:6
        groupVals = valsThisAxis(double(controllerCat(1:n)) == g);
        q1 = prctile(groupVals, 25);
        q3 = prctile(groupVals, 75);
        iqrVal = q3 - q1;
        upperFence = q3 + 1.5*iqrVal;
        whiskerTop = max(groupVals(groupVals <= upperFence));
        
        text(g + xOffset, min(whiskerTop + 0.5, 32), letters{g}, ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold', ...
            'Color', axisColors(a,:))
    end
end

hold off
set(gca, 'FontSize', 16)
ylabel('Interaction force (N)', 'FontSize', 18)
ylim([0, 22])
legend({'X','Y','Z'}, 'Location', 'best')
%% Time length
%groups = [1,2,3,4,5,6];

reaching_time = [];

subj = [];
controller = [];
goal = [];

for j = 1:modes
    for i = 1:movements
        data_all_trials = [];
        for t = 1:trials
            for s = 1:subjects
                reaching_time = [reaching_time; mean(data{j,i,t,s}(end,1))];
                subj = [subj; s];
                controller = [controller; j];
                goal = [goal; i];
         
            end
        end
    end
end

[p, tbl, stats] = anovan(reaching_time, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
[results_controller, m, h, gnames] = multcompare(stats, 'Dimension', 2);

groupMeans = m(:,1);
groupSE = m(:,2);

% 95% CI half-width using the t-distribution (matches the CI convention multcompare uses)
alpha = 0.05;
critVal = tinv(1 - alpha/2, stats.dfe);
ciHalfWidth = critVal * groupSE;

customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};

% Reorder: original index order 1,6,2,3,4,5
newOrder = [1, 6, 2, 3, 4, 5];
groupMeans = groupMeans(newOrder);
ciHalfWidth = ciHalfWidth(newOrder);
customLabels = customLabels(newOrder);

figure
bar(groupMeans)
hold on
errorbar(1:length(groupMeans), groupMeans, ciHalfWidth, 'k', 'LineStyle', 'none', 'LineWidth', 1.2)
hold off



set(gca, 'XTick', 1:length(groupMeans), 'XTickLabel', customLabels,'FontSize', 14)
ylabel('Mean reaching time (s)','FontSize', 16)
ylim([2,3.5])
%title('Mean interaction force by Controller (95% CI)')

% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
figure
results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);


customLabels = {'Force-amplification 1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'Force-amplification 2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

figure
hold on
boxchart(controllerCat, reaching_time, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5])

jitterWidth = 0.3;
xNumeric = double(controllerCat) + (rand(length(controllerCat),1)-0.5)*jitterWidth;
scatter(xNumeric, reaching_time, 4, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.3)

% Compact letter display: FA2 differs from all other groups; the rest share a letter
% Display order: FA1, FA2, ProMP var, Kalman var, ProMP fixed, Kalman fixed
letterLabels = {'A', 'B', 'A', 'A', 'A', 'A'};

yMaxPerGroup = zeros(1,6);
for g = 1:6
    yMaxPerGroup(g) = max(reaching_time(double(controllerCat) == g));
end
letterY = min(yMaxPerGroup + 0.1, 15);   % adjust offset to your data's scale

for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold')
end

hold off
set(gca, 'FontSize', 14)
ylabel('Mean reaching time (s)', 'FontSize', 16)
ylim([0, 16])   % extended slightly above your original [2,3.5] to fit letters

%% Reaching time new graph
customLabels = {'FA-1', 'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed', 'FA-2'};
newOrder = [1, 6, 2, 3, 4, 5];
controllerCat = categorical(controller, [1 6 2 3 4 5], customLabels(newOrder));

% Collapse to one mean per subject per controller
subjControllerPairs = unique([subj, controller], 'rows');
n_pairs = size(subjControllerPairs, 1);
subjectMeans = zeros(n_pairs, 1);
pairControllerCat = categorical(subjControllerPairs(:,2), [1 6 2 3 4 5], customLabels(newOrder));

for k = 1:n_pairs
    thisSubj = subjControllerPairs(k,1);
    thisController = subjControllerPairs(k,2);
    mask = (subj == thisSubj) & (controller == thisController);
    subjectMeans(k) = mean(reaching_time(mask));
end

% Box plot + one jittered point per subject (mean only, no error bars)
figure
hold on
boxchart(controllerCat, reaching_time, 'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5],'LineWidth', 1.5)

jitterWidth = 0.3;
xNumeric = double(pairControllerCat) + (rand(n_pairs,1)-0.5)*jitterWidth;
scatter(xNumeric, subjectMeans, 24, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6)

% Compact letter display, placed just above each group's upper whisker
letterLabels = {'A', 'B', 'A', 'A', 'A', 'A'};

whiskerTop = zeros(1,6);
for g = 1:6
    groupVals = reaching_time(double(controllerCat) == g);
    q1 = prctile(groupVals, 25);
    q3 = prctile(groupVals, 75);
    iqrVal = q3 - q1;
    upperFence = q3 + 1.5*iqrVal;
    whiskerTop(g) = max(groupVals(groupVals <= upperFence));   % actual whisker end, not raw max
end

letterY = min(whiskerTop + 1, 15);
for g = 1:6
    text(g, letterY(g), letterLabels{g}, 'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'FontWeight', 'bold')
end

hold off
set(gca, 'FontSize', 16)
ylabel('Reaching time (s)', 'FontSize', 18)
ylim([0, 7])


%% User data
%2-way ANOVA is not the right tool for the job
perceived_effort(1:modes) = [5,4,4,4,3,5]; %31
comfort(1:modes) = [5,7,6,7,8,6];
perceived_effort(modes+1:2*modes) = [4,4,4,3,3,2]; %32
comfort(modes+1:2*modes) = [8,6,7,7,6,9];
perceived_effort(2*modes+1:3*modes) = [6,5,4,4,3,6]; %33
comfort(2*modes+1:3*modes) = [5,6,6,6,7,4];
perceived_effort(3*modes+1:4*modes) = [4,4,5,5,6,3]; %34
comfort(3*modes+1:4*modes) = [7,5,6,3,4,7];
perceived_effort(4*modes+1:5*modes) = [7,4,4,5,2,7]; %35
comfort(4*modes+1:5*modes) = [3,6,7,5,8,3];
perceived_effort(5*modes+1:6*modes) = [7,8,5,8,5,5]; %36
comfort(5*modes+1:6*modes) = [6,4,7,6,8,6];
perceived_effort(6*modes+1:7*modes) = [8,2,5,3,3,8]; %37
comfort(6*modes+1:7*modes) = [3,8,6,7,7,3];
perceived_effort(7*modes+1:8*modes) = [5,4,4,4,4,4]; %38
comfort(7*modes+1:8*modes) = [5,6,5,6,6,6];
perceived_effort(8*modes+1:9*modes) = [7,3,4,4,3,7]; %39
comfort(8*modes+1:9*modes) = [4,7,6,7,8,5];
perceived_effort(9*modes+1:10*modes) = [6,6,6,6,6,5]; %40
comfort(9*modes+1:10*modes) = [6,6,7,5,6,7];
perceived_effort(10*modes+1:11*modes) = [6,7,6,5,7,5]; %41
comfort(10*modes+1:11*modes) = [7,7,7,8,7,9];
perceived_effort(11*modes+1:12*modes) = [3,4,3,5,4,6]; %42
comfort(11*modes+1:12*modes) = [7,6,7,5,6,6];
perceived_effort(12*modes+1:13*modes) = [6,3,4,3,4,6]; %43
comfort(12*modes+1:13*modes) = [3,8,8,7,6,4];
perceived_effort(13*modes+1:14*modes) = [7,4,6,4,5,6]; %44
comfort(13*modes+1:14*modes) = [4,6,5,6,5,3];
perceived_effort(14*modes+1:15*modes) = [4,4,5,3,5,3]; %45
comfort(14*modes+1:15*modes) = [8,9,7,8,6,8];

subj = [];
controller = [];

for i=1:subjects
    for j=1:modes
        subj = [subj; i];
        controller = [controller; j];
    end
end    

% We need MANOVA and then ANOVA for each one.
%[p, tbl, stats] = anovan(transpose(perceived_effort), {subj, controller},'model', 'interaction','varnames', {'Subject','Controller'})

[d, p, stats] = manova1([transpose(perceived_effort), transpose(comfort)], controller)

Y = [transpose(perceived_effort), transpose(comfort)];
groupMeans = grpstats(Y, controller);
disp(groupMeans)

%%

customLabels = {'Admittance', 'ProMP var', 'Kalman var', 'ProMP fixed', 'Kalman fixed', 'Admittance 2'}; % adjust to your 6 labels
n_controllers = 6;

% X-position for each controller's effort column and comfort column
% Effort at positions 1,3,5,7,9,11 ; Comfort at positions 2,4,6,8,10,12
xEffort = 2*controller - 1;
xComfort = 2*controller;

% Add small horizontal jitter so overlapping subject points are visible
jitterAmount = 0.12;
xEffort_jit = xEffort + (rand(size(xEffort))-0.5)*jitterAmount;
xComfort_jit = xComfort + (rand(size(xComfort))-0.5)*jitterAmount;

figure
hold on
scatter(xEffort_jit, transpose(perceived_effort), 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6)
scatter(xComfort_jit, transpose(comfort), 50, 'r', 'filled', 'MarkerFaceAlpha', 0.6)
hold off

% Build x-tick labels: one pair of ticks per controller
xticks(1:2*n_controllers)
xTickLabels = cell(1, 2*n_controllers);
for c = 1:n_controllers
    xTickLabels{2*c-1} = sprintf('%s\nEffort', customLabels{c});
    xTickLabels{2*c}   = sprintf('%s\nComfort', customLabels{c});
end
xticklabels(xTickLabels)
xtickangle(45)

ylabel('Rating')
title('Perceived Effort and Comfort by Controller (individual subjects)')
legend({'Perceived Effort', 'Comfort'}, 'Location', 'best')
xlim([0.5, 2*n_controllers+0.5])


%% Current
%groups = [1,2,3,4,5,6];
Current_mean = [];

subj = [];
controller = [];
goal = [];

for j = 1:modes
    for i = 1:movements
        data_all_trials = [];
        for t = 1:trials
            for s = 1:subjects
                Current_mean = [Current_mean; mean(sqrt(data{j,i,t,s}(1:end-1,37).^2+...
                    data{j,i,t,s}(1:end-1,38).^2+ data{j,i,t,s}(1:end-1,39).^2 + data{j,i,t,s}(1:end-1,40).^2 +...
                    data{j,i,t,s}(1:end-1,41).^2+ data{j,i,t,s}(1:end-1,42).^2+ data{j,i,t,s}(1:end-1,43).^2))];
                subj = [subj; s];
                controller = [controller; j];
                goal = [goal; i];
         
            end
        end
    end
end

[p, tbl, stats] = anovan(Current_mean, {subj, controller, goal},'model', 'interaction','varnames', {'Subject','Controller','Target'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for Goal (dimension 3)
figure
results_subject = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);


%%%
%% Build mean prediction error per subject/controller/target
predError_mean = [];
predError_meanX = [];
predError_meanY = [];
predError_meanZ = [];
subj = [];
controller = [];
goal = [];

prediction_modes = [2,3,4,5];

for j = prediction_modes
    for i = 1:movements
        for t = 1:trials
            for s = 1:subjects
                current_data = data{j,i,t,s};
                time = current_data(:,1);
                
                predErrorX = zeros(length(time),1);
                predErrorY = zeros(length(time),1);
                predErrorZ = zeros(length(time),1);
                
                for k = 1:length(time)
                    predX = current_data(k,20);
                    predY = current_data(k,21);
                    predZ = current_data(k,22);
                    s_diff = current_data(k,32);
                    
                    [~, idx] = min(abs(time - (time(k)+s_diff)));
                    
                    truePosX = current_data(idx,2);
                    truePosY = current_data(idx,3);
                    truePosZ = current_data(idx,4);
                    
                    predErrorX(k) = abs(predX-truePosX);
                    predErrorY(k) = abs(predY-truePosY);
                    predErrorZ(k) = abs(predZ-truePosZ);
                end
                
                % Total (Euclidean) prediction error per sample, then mean over the trial
                totalError = sqrt(predErrorX.^2 + predErrorY.^2 + predErrorZ.^2);
                trial_mean_error = mean(totalError(2:end), 'omitnan');
                
                predError_mean = [predError_mean; trial_mean_error];
                predError_meanX = [predError_meanX; mean(predErrorX)];
                predError_meanY = [predError_meanY; mean(predErrorY)];
                predError_meanZ = [predError_meanZ; mean(predErrorZ)];
                subj = [subj; s];
                controller = [controller; j];
                goal = [goal; i];
            end
        end
    end
end




%% Two-way (three-factor) ANOVA
[p, tbl, stats] = anovan(predError_mean, {subj, controller, goal}, ...
    'model', 'interaction', ...
    'varnames', {'Subject','Controller','Target'});

% Pairwise comparisons for Controller (dimension 2)
figure
results_controller = multcompare(stats, 'Dimension', 2)
% Pairwise comparisons for Target/Goal (dimension 3)



figure
results_goal = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% Bar chart of mean prediction error by controller, with 95% CI
controllerLevels = unique(controller);
customLabels = {'ProMP var', 'UKF var', 'ProMP fixed', 'UKF fixed'}; % adjust order to match modes [2,3,4,5]

alpha = 0.05;
means = zeros(length(controllerLevels), 1);
ciHalfWidth = zeros(length(controllerLevels), 1);

for c = 1:length(controllerLevels)
    vals = predError_mean(controller == controllerLevels(c));
    n = length(vals);
    means(c) = mean(vals);
    se = std(vals) / sqrt(n);
    tCrit = tinv(1 - alpha/2, n-1);
    ciHalfWidth(c) = tCrit * se;
end

figure
bar(means, 'FaceColor', [0 0.2 0.55], 'FaceAlpha', 0.6)
hold on
%errorbar(1:length(means), means, ciHalfWidth, 'Color', [0 0.2 0.55], 'LineStyle', 'none', 'LineWidth', 1.5)
errorbar(1:length(means), means, ciHalfWidth, 'k', 'LineStyle', 'none', 'LineWidth', 1.5)
% Significance brackets: ProMP variable vs UKF variable, and ProMP variable vs UKF fixed
yMaxAll = max(means + ciHalfWidth);
yRange = yMaxAll - min(means - ciHalfWidth);

% ProMP variable (1) vs UKF variable (2)
x1 = 1; x2 = 2;
yBracket1 = yMaxAll + 0.08*yRange;
plot([x1 x1 x2 x2], [yBracket1-0.02*yRange, yBracket1, yBracket1, yBracket1-0.02*yRange], 'k-', 'LineWidth', 1.5)
text(mean([x1 x2]), yBracket1 + 0.02*yRange, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')

% ProMP variable (1) vs UKF fixed (4)
x1 = 1; x2 = 4;
yBracket2 = yMaxAll + 0.20*yRange;
plot([x1 x1 x2 x2], [yBracket2-0.02*yRange, yBracket2, yBracket2, yBracket2-0.02*yRange], 'k-', 'LineWidth', 1.5)
text(mean([x1 x2]), yBracket2 + 0.02*yRange, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')

hold off
set(gca, 'XTick', 1:length(means), 'XTickLabel', customLabels,'FontSize', 18)
ylabel('Error (m)','FontSize', 20)
ylim([min(means - ciHalfWidth) - 0.05*yRange, yMaxAll + 0.30*yRange])
%title('Mean Prediction Error by Controller (95% CI)')

%%
% Find which specific (j,i,t,s) combinations produced NaN or empty data
for j = prediction_modes
    for i = 1:movements
        for t = 1:trials
            for s = 1:subjects
                current_data = data{j,i,t,s};
                if isempty(current_data)
                    continue
                end
                
                cols_to_check = [2,3,4,20,21,22,32];
                nanCount = sum(any(isnan(current_data(:, cols_to_check)), 2));
                
                if nanCount > 0
                    fprintf('mode=%d, movement=%d, trial=%d, subject=%d: %d rows with NaN\n', ...
                        j, i, t, s, nanCount);
                end
            end
        end
    end
end


%%
[p, tbl, stats] = anovan(predError_meanX, {subj, controller, goal}, ...
    'model', 'interaction', ...
    'varnames', {'Subject','Controller','Target'});

% Pairwise comparisons for Controller (dimension 2)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Target/Goal (dimension 3)
figure
results_goal = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);
%%
[p, tbl, stats] = anovan(predError_meanY, {subj, controller, goal}, ...
    'model', 'interaction', ...
    'varnames', {'Subject','Controller','Target'});

% Pairwise comparisons for Controller (dimension 2)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Target/Goal (dimension 3)
figure
results_goal = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);
%%
[p, tbl, stats] = anovan(predError_meanZ, {subj, controller, goal}, ...
    'model', 'interaction', ...
    'varnames', {'Subject','Controller','Target'});

% Pairwise comparisons for Controller (dimension 2)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Target/Goal (dimension 3)
figure
results_goal = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);



%%
%data{1,1} = load("new_table/subject18mode3mov1.csv");
%data{2,1} = load("var_horizon1/r2.csv");
%data{3,1} = load("var_horizon1/r3.csv");

i = 1;
s_diff = 0.030;



%%
% 
% figure()
% subplot(3,1,1)
% plot(data{i,1}(:,1),data{i,1}(:,2))
% title('posX')
% 
% subplot(3,1,2)
% plot(data{i,1}(:,1),data{i,1}(:,3))
% title('posY')
% 
% subplot(3,1,3)
% plot(data{i,1}(:,1),data{i,1}(:,4))
% title('posZ')
% 
% figure()
% subplot(3,1,1)
% plot(data{i,1}(:,1),data{i,1}(:,8))
% title('VelX')
% 
% subplot(3,1,2)
% plot(data{i,1}(:,1),data{i,1}(:,9))
% title('VelY')
% 
% subplot(3,1,3)
% plot(data{i,1}(:,1),data{i,1}(:,10))
% title('VelZ')
% 
% figure
% subplot(3,1,1)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,8)./max(data{i,1}(:,8)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio VelX')
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,9)./max(data{i,1}(:,9)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio VelY')
% hold off
% 
% subplot(3,1,3)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,10)./max(abs(data{i,1}(:,10))))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(abs(data{i,1}(:,36))))
% title('ratio VelZ')
% hold off
% 
% % figure()
% % hold on
% % plot(data{i,1}(:,1),data{i,1}(:,2))
% % plot(data{i,1}(:,1),data{i,1}(:,3))
% % plot(data{i,1}(:,1),data{i,1}(:,4))
% % plot(data{i,1}(:,1),data{i,1}(:,5))
% % plot(data{i,1}(:,1),data{i,1}(:,6))
% % plot(data{i,1}(:,1),data{i,1}(:,7))
% % hold off
% % legend("posX","posY","posZ","rot1","rot2","rot3")
% % title('True trajectory')
% % 
% % 
% % figure()
% % hold on
% % plot(data{i,1}(:,1),data{i,1}(:,8))
% % plot(data{i,1}(:,1),data{i,1}(:,9))
% % plot(data{i,1}(:,1),data{i,1}(:,10))
% % plot(data{i,1}(:,1),data{i,1}(:,11))
% % plot(data{i,1}(:,1),data{i,1}(:,12))
% % plot(data{i,1}(:,1),data{i,1}(:,13))
% % hold off
% % legend("velX","velY","velZ","velrot1","velrot2","velrot3")
% % title('True velocity')
% 
% figure()
% plot(data{i,1}(:,1), sqrt((data{i,1}(:,8).^2+data{i,1}(:,9).^2+data{i,1}(:,10).^2)))
% title('Total velocity')
% 
% 
% figure()
% subplot(3,1,1)
% plot(data{i,1}(:,1),data{i,1}(:,14))
% title('TorqueX')
% 
% subplot(3,1,2)
% plot(data{i,1}(:,1),data{i,1}(:,15))
% title('TorqueY')
% 
% subplot(3,1,3)
% plot(data{i,1}(:,1),data{i,1}(:,16))
% title('TorqueZ')
% 
% figure
% subplot(3,1,1)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,14)./max(data{i,1}(:,14)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueX')
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,15)./max(data{i,1}(:,15)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueY')
% hold off
% 
% 
% subplot(3,1,3)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,16)./max(abs(data{i,1}(:,16))))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(abs(data{i,1}(:,36))))
% title('ratio TorqueZ')
% hold off
% 
% 
% figure()
% subplot(3,1,1)
% 
% plot(data{i,1}(:,1),data{i,1}(:,14).*data{i,1}(:,8))
% title('TorqueX*VelocityX')
% 
% 
% subplot(3,1,2)
% plot(data{i,1}(:,1),data{i,1}(:,15).*data{i,1}(:,9))
% title('TorqueY*VelocityY')
% 
% subplot(3,1,3)
% plot(data{i,1}(:,1),data{i,1}(:,16).*data{i,1}(:,10))
% title('TorqueZ*VelocityZ')
% 
% figure
% subplot(3,1,1)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,14).*data{i,1}(:,8)./max(data{i,1}(:,14).*data{i,1}(:,8)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueX*VelocityX')
% hold off
% 
% subplot(3,1,2)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,15).*data{i,1}(:,9)./max(data{i,1}(:,15).*data{i,1}(:,9)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueY*VelocityY')
% hold off
% 
% 
% subplot(3,1,3)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,16).*data{i,1}(:,10)./max(data{i,1}(:,16).*data{i,1}(:,10)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueZ*VelocityZ')
% hold off
% 
% %%
% 

%% Calculate the error in prediction position
time = data{i,1}(:,1);
predErrorX = zeros(length(time),1);
predErrorY = zeros(length(time),1);
predErrorZ = zeros(length(time),1);
for j = 1:length(time)
    predX = data{i,1}(j,20);
    predY = data{i,1}(j,21);
    predZ = data{i,1}(j,22);
    s_diff = data{i,1}(j,32);

    [~, idx] = min(abs(time - (time(j)+s_diff)));
    truePosX = data{i,1}(idx,2);
    truePosY = data{i,1}(idx,3);
    truePosZ = data{i,1}(idx,4);
    predErrorX(j) = abs(predX-truePosX);
    predErrorY(j) = abs(predY-truePosY); 
    predErrorZ(j) = abs(predZ-truePosZ); 
end

figure()
subplot(3,1,1)
plot(time(2:end),predErrorX(2:end)*100)
title('Pred error for X in cm')

subplot(3,1,2)
plot(time(2:end),predErrorY(2:end)*100)
title('Pred error for Y in cm')

subplot(3,1,3)
plot(time(2:end),predErrorZ(2:end)*100)
title('Pred error for Z in cm')

%% Calculate the error in prediction Velocity
time = data{i,1}(:,1);
predErrorX = zeros(length(time),1);
predErrorY = zeros(length(time),1);
predErrorZ = zeros(length(time),1);
for j = 1:length(time)
    predX = data{i,1}(j,26);
    predY = data{i,1}(j,27);
    predZ = data{i,1}(j,28);
    s_diff = data{i,1}(j,32);

    [~, idx] = min(abs(time - (time(j)+s_diff)));
    trueVelX = data{i,1}(idx,8);
    trueVelY = data{i,1}(idx,9);
    trueVelZ = data{i,1}(idx,10);
    predErrorX(j) = abs(predX-trueVelX);
    predErrorY(j) = abs(predY-trueVelY); 
    predErrorZ(j) = abs(predZ-trueVelZ); 
end

figure()
subplot(3,1,1)
plot(time(2:end),predErrorX(2:end)*100)
title('Pred velocity error for X in cm')

subplot(3,1,2)
plot(time(2:end),predErrorY(2:end)*100)
title('Pred velocity for Y in cm')

subplot(3,1,3)
plot(time(2:end),predErrorZ(2:end)*100)
title('Pred velocity for Z in cm')

% %
% figure()
% subplot(3,1,1)
% plot(data{i,1}(:,1),data{i,1}(:,33))
% title('Kx')
% 
% subplot(3,1,2)
% plot(data{i,1}(:,1),data{i,1}(:,34))
% title('Ky')
% 
% subplot(3,1,3)
% plot(data{i,1}(:,1),data{i,1}(:,35))
% title('Kz')
% 
% %% sensor reading
% figure()
% plot(data{i,1}(:,1),data{i,1}(:,36))
% title('Human effort')
% 
% 
% %%
% idx = abs(data{i,1}(:,8)) < 0.1;
% figure
% subplot(3,1,1)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,14).*(idx)./max(data{i,1}(:,14)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueX*(1/VelocityX)')
% hold off
% 
% idx = abs(data{i,1}(:,9)) < 0.1;
% subplot(3,1,2)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,15).*idx./max(data{i,1}(:,15)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueY*1/VelocityY')
% hold off
% 
% idx = abs(data{i,1}(:,10)) < 0.1;
% subplot(3,1,3)
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,16).*idx./max(data{i,1}(:,16)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% title('ratio TorqueZ*1/VelocityZ')
% hold off
% 
% %%
% figure()
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,8)./max(data{i,1}(:,8)))
% plot(data{i,1}(:,1),data{i,1}(:,14)./max(data{i,1}(:,14)))
% plot(data{i,1}(:,1),data{i,1}(:,33)./max(data{i,1}(:,33)))
% plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
% legend('vel', 'force','stiffness','effort')
% hold off
% 
% figure()
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,2))
% plot(data{i,1}(:,1),data{i,1}(:,20))
% plot(data{i,1}(:,1),sign(data{i,1}(:,8)))
% plot(data{i,1}(:,1),data{i,1}(:,8))
% 
% hold off