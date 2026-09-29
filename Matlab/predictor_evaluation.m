modes = [1,6];
movements = 5;
trials = 5;
subjects = 15;

data_ProMP = cell(length(modes),trials,movements,subjects);
data_kalman = cell(length(modes),trials,movements,subjects);
data = cell(length(modes),trials,movements,subjects);
for s = 1:subjects
    for j = modes
        for i = 1:trials
            for z = 1:movements

                data{j,i,z,s} = readmatrix(sprintf("predictor_validation/subject%dmode%dmov%d%d.csv", s+30, j,z,i));
                data_ProMP{j,i,z,s} = readmatrix(sprintf("predictor_validation/ProMP_subject%dmode%dmov%d%d.csv", s+30, j,z,i));
                data_kalman{j,i,z,s} = readmatrix(sprintf("predictor_validation/Kalman_subject%dmode%dmov%d%d.csv", s+30, j,z,i));

                fprintf('%d %d %d %d\n', s, j, i, z)
            end
        end
    end
end


%%
axisCols_base = [2, 3, 4];   % columns in `data{...}` for X, Y, Z ground truth
predictorAxisCols = {[3, 5, 7], [3, 4, 5]};   % per-predictor pred-column sets
axisNames = {'X','Y','Z'};

predictorsData = {data_ProMP, data_kalman};
predictorNames = {'ProMP','Kalman'};

% Preallocate and loop
maxIters = numel(predictorsData) * length(modes) * trials * movements * subjects * 3;
resultsCell = cell(maxIters, 1);
counter = 0;

for p = 1:numel(predictorsData)
    predData = predictorsData{p};
    pname = predictorNames{p};
    axisCols_pred = predictorAxisCols{p};   
    
    for jIdx = 1:length(modes)
        j = modes(jIdx);
        for s = 1:subjects
            for i = 1:trials
                for z = 1:movements
                
                    fprintf('%d %d %d %d\n', j, i,z,s)
                    baseData = data{j,i,z,s};
                    predAll = predData{j,i,z,s};
                    
                    if isempty(baseData) || isempty(predAll)
                        continue
                    end
                    
                    for axisIdx = 1:3
                        A = predAll(:, [1, 2, axisCols_pred(axisIdx)]);
                        baseAxis = baseData(:, [1, axisCols_base(axisIdx)]);
                        
                        [horizon, errorVal, errorInt] = computeAxisMetrics(A, baseAxis);
                        nLabels = size(horizon, 2);
                        if nLabels == 0
                            continue
                        end
                        
                        counter = counter + 1;
                        resultsCell{counter} = table( ...
                            s, j, z, i, ...
                            string(pname), string(axisNames{axisIdx}), ...
                            mean(horizon(1,:)), mean(horizon(2,:)), mean(horizon(3,:)), ...
                            mean(errorVal(1,:)), mean(errorVal(2,:)), mean(errorVal(3,:)), ...
                            mean(errorInt(1,:)), mean(errorInt(2,:)), mean(errorInt(3,:)), ...
                            'VariableNames', {'Subject','Controller','Movement','Trial','Predictor','Axis', ...
                            'Horizon10cm','Horizon5cm','Horizon1cm', ...
                            'Error500ms','Error100ms','Error50ms', ...
                            'ErrorInt500ms','ErrorInt100ms','ErrorInt50ms'});
                    end
                end
            end
        end
    end
end

resultsCell = resultsCell(1:counter);
resultsTable = vertcat(resultsCell{:});

% Local function — computes horizon and error metrics for one axis
function [horizon, errorVal, errorInt] = computeAxisMetrics(A, baseAxisData)
    labels = unique(A(:,1));
    nLabels = length(labels);
    horizon = zeros(3, nLabels);
    errorVal = zeros(3, nLabels);
    errorInt = zeros(3, nLabels);
    
    for k = 1:nLabels
        grp = A(A(:,1) == labels(k), :);

        [~, uniqueIdx] = unique(grp(:,2), 'stable');
        grp = grp(uniqueIdx, :);

        baseCut = baseAxisData(baseAxisData(:,1) >= grp(1,2), :);

        if size(grp,1) <= 1 || isempty(baseCut)
            continue
        end

        predInterp = spline(grp(:,2), grp(:,3), baseCut(:,1));
        if length(predInterp) <= 1
            continue
        end
        
        t0 = baseCut(1,1);
        tEndPred = grp(end,2);
        errAll = abs(predInterp(:) - baseCut(:,2));
        
        % --- Horizon: cap at samples within the prediction's valid time range ---
        validMask = baseCut(:,1) <= tEndPred;
        lastValidIdx = find(validMask, 1, 'last');
        if isempty(lastValidIdx)
            lastValidIdx = 0;
        end
        
        over10 = find(errAll(1:lastValidIdx) >= 0.10, 1, 'first');
        if isempty(over10)
            runLen = lastValidIdx;
        else
            runLen = over10 - 1;
        end
        
        if runLen >= 1
            horizon(1,k) = baseCut(runLen,1) - t0;
            idx05 = find(errAll(1:runLen) < 0.05, 1, 'last');
            if ~isempty(idx05)
                horizon(2,k) = baseCut(idx05,1) - t0;
                idx01 = find(errAll(1:idx05) < 0.01, 1, 'last');
                if ~isempty(idx01)
                    horizon(3,k) = baseCut(idx01,1) - t0;
                end
            end
        end
        
        % --- Error at fixed time windows (500ms, 100ms, 50ms) ---
        relTime = baseCut(:,1) - t0;
        
        idx500 = find(relTime <= 0.5, 1, 'last');
        if ~isempty(idx500)
            errorVal(1,k) = errAll(idx500);
            errorInt(1,k) = sum(errAll(1:idx500));
        end
        
        idx100 = find(relTime <= 0.1, 1, 'last');
        if ~isempty(idx100)
            errorVal(2,k) = errAll(idx100);
            errorInt(2,k) = sum(errAll(1:idx100));
        end
        
        idx50 = find(relTime <= 0.05, 1, 'last');
        if ~isempty(idx50)
            errorVal(3,k) = errAll(idx50);
            errorInt(3,k) = sum(errAll(1:idx50));
        end
    end
end

save('prediction_results.mat', 'resultsTable')
%% MANOVA
load('prediction_results.mat');  % restores resultsTable to the workspace

X = resultsTable{:, {'Horizon10cm','Horizon5cm', 'Horizon1cm', 'Error500ms', 'Error100ms', 'Error50ms', 'ErrorInt500ms', 'ErrorInt100ms', 'ErrorInt50ms'}};
group = categorical(resultsTable.Predictor);

[d, p, stats] = manova1(X, group)

% ANOVA for each variable
metrics = {'Horizon10cm','Horizon5cm','Horizon1cm','Error500ms','Error100ms','Error50ms','ErrorInt500ms','ErrorInt100ms','ErrorInt50ms'};

pValuesOneWay = zeros(length(metrics), 1);
for m = 1:length(metrics)
    [p, tbl, stats] = anova1(resultsTable.(metrics{m}), resultsTable.Predictor, 'off');
    pValuesOneWay(m) = p;
end

alpha = 0.05;
alphaCorrected = alpha / length(metrics);
significant = pValuesOneWay < alphaCorrected;

table(metrics', pValuesOneWay, significant, 'VariableNames', {'Metric','pValue','SignificantBonferroni'})

%% paired sample t-test for each variable
metrics = {'Horizon10cm','Horizon5cm','Horizon1cm','Error500ms','Error100ms','Error50ms','ErrorInt500ms','ErrorInt100ms','ErrorInt50ms'};

pValuesPaired = zeros(length(metrics), 1);
tStats = zeros(length(metrics), 1);
nPairs = zeros(length(metrics), 1);

keyVars = {'Subject','Controller','Movement','Trial','Axis'};

for m = 1:length(metrics)
    subTable = resultsTable(:, [keyVars, {'Predictor', metrics{m}}]);
    wideTable = unstack(subTable, metrics{m}, 'Predictor');
    
    % Drop rows where either predictor is missing (unmatched trials)
    wideTable = rmmissing(wideTable);
    
    [~, p, ~, stats] = ttest(wideTable.Kalman, wideTable.ProMP);
    
    pValuesPaired(m) = p;
    tStats(m) = stats.tstat;
    nPairs(m) = height(wideTable);
end

alpha = 0.05;
alphaCorrected = alpha / length(metrics);
significant = pValuesPaired < alphaCorrected;

table(metrics', pValuesPaired, tStats, nPairs, significant, ...
    'VariableNames', {'Metric','pValue','tStat','N_pairs','SignificantBonferroni'})
%% Plot of horizons
predictors = {'Kalman', 'ProMP'};
metrics = {'Horizon10cm', 'Horizon5cm', 'Horizon1cm'};

alpha = 0.05;
means = zeros(length(predictors), length(metrics));
ciHalfWidth = zeros(length(predictors), length(metrics));

for pIdx = 1:length(predictors)
    for mIdx = 1:length(metrics)
        vals = resultsTable.(metrics{mIdx})(resultsTable.Predictor == predictors{pIdx});
        n = length(vals);
        
        means(pIdx, mIdx) = mean(vals);
        se = std(vals) / sqrt(n);
        tCrit = tinv(1 - alpha/2, n-1);
        ciHalfWidth(pIdx, mIdx) = tCrit * se;
    end
end

figure
b = bar(means);
hold on

for mIdx = 1:length(metrics)
    xPos = b(mIdx).XEndPoints;
    errorbar(xPos, means(:,mIdx), ciHalfWidth(:,mIdx), 'k', 'LineStyle', 'none', 'LineWidth', 1.2)
end

predictors = {'UKF', 'ProMP'};
set(gca, 'XTickLabel', predictors,'FontSize', 14)
ylabel('Time (s)','FontSize', 16)
metrics = {'Horizon 10cm', 'Horizon 5cm', 'Horizon 1cm'};
%title('Prediction Horizon by Predictor (95% CI)')
legend(metrics, 'Location', 'best','FontSize', 14)

%% Plot of horizons
predictorLabels = {'Kalman', 'ProMP'};
metrics = {'Horizon10cm', 'Horizon5cm', 'Horizon1cm'};
alpha = 0.05;
means = zeros(length(predictorLabels), length(metrics));
ciHalfWidth = zeros(length(predictorLabels), length(metrics));
for pIdx = 1:length(predictorLabels)
    for mIdx = 1:length(metrics)
        vals = resultsTable.(metrics{mIdx})(resultsTable.Predictor == predictorLabels{pIdx});
        n = length(vals);
        means(pIdx, mIdx) = mean(vals);
        se = std(vals) / sqrt(n);
        tCrit = tinv(1 - alpha/2, n-1);
        ciHalfWidth(pIdx, mIdx) = tCrit * se;
    end
end

%[0 0 0; 0.4 0.7 0.9; 0.9 0.4 0.3]

figure
%b = bar(means);
b = bar(means, 'FaceAlpha', 0.6);
b(1).FaceColor = [0 0.2 0.55];       % blue
b(2).FaceColor = [0.75 0.35 0.05];   % darker orange
b(3).FaceColor = [0.75 0.6 0.05];    % darker yellow
hold on
for mIdx = 1:length(metrics)
    xPos = b(mIdx).XEndPoints;
    errorbar(xPos, means(:,mIdx), ciHalfWidth(:,mIdx), 'k', 'LineStyle', 'none', 'LineWidth', 1.5)
end

% Significance brackets: Kalman vs ProMP, for each of the three metrics
yMaxAll = max(means(:) + ciHalfWidth(:));
yRange = yMaxAll - min(means(:) - ciHalfWidth(:));

for mIdx = 1:length(metrics)
    x1 = b(mIdx).XEndPoints(1);   % Kalman bar for this metric
    x2 = b(mIdx).XEndPoints(2);   % ProMP bar for this metric
    
    yBracket = yMaxAll + 0.0*yRange + (length(metrics)-mIdx)*0.10*yRange;
    
    plot([x1 x1 x2 x2], [yBracket-0.02*yRange, yBracket, yBracket, yBracket-0.02*yRange], 'k-', 'LineWidth', 1.5)
    text(mean([x1 x2]), yBracket + 0.02*yRange, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')
end

xTickLabelsFinal = {'UKF', 'ProMP'};
set(gca, 'XTickLabel', xTickLabelsFinal,'FontSize', 20)
ylabel('Time (s)','FontSize', 20)
metricLabels = {'Horizon 10 cm', 'Horizon 5 cm', 'Horizon 1 cm'};
legend(metricLabels, 'Location', 'best','FontSize', 20)
ylim([0, yMaxAll + 0.40*yRange])

hold off

%% Plot of errors (im not convinced by the logarithmic scale)
predictors = {'Kalman', 'ProMP'};
metrics = {'Error500ms', 'Error100ms', 'Error50ms'};
alpha = 0.05;
means = zeros(length(predictors), length(metrics));
ciHalfWidth = zeros(length(predictors), length(metrics));
for pIdx = 1:length(predictors)
    for mIdx = 1:length(metrics)
        vals = resultsTable.(metrics{mIdx})(resultsTable.Predictor == predictors{pIdx});
        n = length(vals);
        means(pIdx, mIdx) = mean(vals, 'omitnan');
        se = std(vals, 'omitnan') / sqrt(n);
        tCrit = tinv(1 - alpha/2, n-1);
        ciHalfWidth(pIdx, mIdx) = tCrit * se;
    end
end

figure
b = bar(means, 'FaceAlpha', 0.6);
b(1).FaceColor = [0 0.2 0.55];       % blue
b(2).FaceColor = [0.75 0.35 0.05];   % darker orange
b(3).FaceColor = [0.75 0.6 0.05];    % darker yellow
hold on
for mIdx = 1:length(metrics)
    xPos = b(mIdx).XEndPoints;
    errorbar(xPos, means(:,mIdx), ciHalfWidth(:,mIdx), 'k', 'LineStyle', 'none', 'LineWidth', 1.5)
end


% Significance brackets: Kalman vs ProMP, for each of the three metrics
yMaxAll = max(means(:) + ciHalfWidth(:));
yRange = yMaxAll - min(means(:) - ciHalfWidth(:));

for mIdx = 1:length(metrics)
    x1 = b(mIdx).XEndPoints(1);   % Kalman bar for this metric
    x2 = b(mIdx).XEndPoints(2);   % ProMP bar for this metric
    
    yBracket = yMaxAll + 0.2*yRange + (length(metrics)-mIdx)*0.20*yRange;
    
    plot([x1 x1 x2 x2], [yBracket-0.05*yRange, yBracket, yBracket, yBracket-0.05*yRange], 'k-', 'LineWidth', 1.5)
    text(mean([x1 x2]), yBracket + 0.02*yRange, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')
end

predictorLabels = {'UKF', 'ProMP'};
set(gca, 'XTickLabel', predictorLabels)
set(gca, 'YScale', 'log','FontSize', 20)
ylabel('Error (m)','FontSize', 20)
metricLabels = {'Error 500 ms', 'Error 100 ms', 'Error 50 ms'};
legend(metricLabels, 'Location', 'best','FontSize', 20)
hold off

%% Plot of integrated errors (im not convinced by the logarithmic scale)
predictors = {'Kalman', 'ProMP'};
metrics = {'ErrorInt500ms', 'ErrorInt100ms', 'ErrorInt50ms'};
alpha = 0.05;
means = zeros(length(predictors), length(metrics));
ciHalfWidth = zeros(length(predictors), length(metrics));
for pIdx = 1:length(predictors)
    for mIdx = 1:length(metrics)
        vals = resultsTable.(metrics{mIdx})(resultsTable.Predictor == predictors{pIdx});
        n = length(vals);
        means(pIdx, mIdx) = mean(vals, 'omitnan');
        
        se = std(vals, 'omitnan') / sqrt(n);
        tCrit = tinv(1 - alpha/2, n-1);
        ciHalfWidth(pIdx, mIdx) = tCrit * se;
    end
end

figure
b = bar(means, 'FaceAlpha', 0.6);
b(1).FaceColor = [0 0.2 0.55];       % blue
b(2).FaceColor = [0.75 0.35 0.05];   % darker orange
b(3).FaceColor = [0.75 0.6 0.05];    % darker yellow
hold on
for mIdx = 1:length(metrics)
    xPos = b(mIdx).XEndPoints;
    errorbar(xPos, means(:,mIdx), ciHalfWidth(:,mIdx), 'k', 'LineStyle', 'none', 'LineWidth', 1.5)
end

% Significance brackets: Kalman vs ProMP, for each of the three metrics
yMaxAll = max(means(:) + ciHalfWidth(:));
yRange = yMaxAll - min(means(:) - ciHalfWidth(:));

for mIdx = 1:length(metrics)
    x1 = b(mIdx).XEndPoints(1);   % Kalman bar for this metric
    x2 = b(mIdx).XEndPoints(2);   % ProMP bar for this metric
    
    yBracket = yMaxAll + 0.2*yRange + (length(metrics)-mIdx)*0.5*yRange;
    
    plot([x1 x1 x2 x2], [yBracket-0.10*yRange, yBracket, yBracket, yBracket-0.10*yRange], 'k-', 'LineWidth', 1.5)
    text(mean([x1 x2]), yBracket + 0.02*yRange, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold')
end

predictorLabels = {'UKF', 'ProMP'};
set(gca, 'XTickLabel', predictorLabels)
set(gca, 'YScale', 'log','FontSize', 20)
ylabel('Integrated error (m)','FontSize', 20)
metricLabels = {'Int. error 500 ms', 'Int. error 100 ms', 'Int. error 50 ms'};
legend(metricLabels, 'Location', 'best','FontSize', 18)
hold off



%% 
[p, tbl, stats] = anovan(resultsTable{:, {'Horizon10cm'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}},'model', 'interaction','varnames', {'Subject','Predictor'})
% Pairwise comparisons for Predictor (dimension 2 in the anovan call)
figure
results_predictor_10cm = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
% figure
% results_subject = multcompare(stats, 'Dimension', 1);
% figure
% results_axis = multcompare(stats, 'Dimension', 3);
% % Pairwise comparisons for each Subject-Controller combination (interaction)
% figure
% results_interaction = multcompare(stats, 'Dimension', [1 2]);

%
[p, tbl, stats] = anovan(resultsTable{:, {'Horizon5cm'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}},'model', 'interaction','varnames', {'Subject','Predictor'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_predictor_5cm = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
% figure
% results_subject = multcompare(stats, 'Dimension', 1);
% figure
% results_axis = multcompare(stats, 'Dimension', 3);
% % Pairwise comparisons for each Subject-Controller combination (interaction)
% figure
% results_interaction = multcompare(stats, 'Dimension', [1 2]);

%
[p, tbl, stats] = anovan(resultsTable{:, {'Horizon1cm'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}},'model', 'interaction','varnames', {'Subject','Predictor'})
% Pairwise comparisons for Predictor (dimension 2 in the anovan call)
figure
results_predictor_1cm = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
% figure
% results_subject = multcompare(stats, 'Dimension', 1);
% figure
% results_axis = multcompare(stats, 'Dimension', 3);
% % Pairwise comparisons for each Subject-Controller combination (interaction)
% figure
% results_interaction = multcompare(stats, 'Dimension', [1 2]);






%% 
[p, tbl, stats] = anovan(resultsTable{:, {'Error500ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% 
[p, tbl, stats] = anovan(resultsTable{:, {'Error100ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% 
[p, tbl, stats] = anovan(resultsTable{:, {'Error50ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% 
[p, tbl, stats] = anovan(resultsTable{:, {'ErrorInt500ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% 
[p, tbl, stats] = anovan(resultsTable{:, {'ErrorInt100ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

%% 
[p, tbl, stats] = anovan(resultsTable{:, {'ErrorInt50ms'}}, {resultsTable{:, {'Subject'}}, resultsTable{:, {'Predictor'}}, resultsTable{:, {'Axis'}}},'model', 'interaction','varnames', {'Subject','Predictor', 'Axis'})
% Pairwise comparisons for Controller (dimension 2 in the anovan call)
figure
results_controller = multcompare(stats, 'Dimension', 2);
% Pairwise comparisons for Subject (dimension 1)
figure
results_subject = multcompare(stats, 'Dimension', 1);
figure
results_axis = multcompare(stats, 'Dimension', 3);
% Pairwise comparisons for each Subject-Controller combination (interaction)
figure
results_interaction = multcompare(stats, 'Dimension', [1 2]);

