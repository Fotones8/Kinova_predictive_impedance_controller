% data {mode, movement}
modes = 6;
movements = 5;
trials = 5;
subject = 45;

data = cell(5,5,5);
for j = 1:modes
    for i = 1:trials
        for z = 1:movements
            data{j,i,z} = readmatrix(sprintf("subject%d/subject%dmode%dmov%d%d.csv", subject, subject, j,z,i));
        %fprintf('%d %d\n', j, i)
        end
    end
end

% %% Clip data based on total velocity threshold
% data_clipped = cell(modes, trials, movements);
% 
% for j = 1:modes
%     for i = 1:trials
%         for z = 1:movements
%             current_data = data{j,i,z};
% 
%             % Compute total velocity magnitude
%             total_velocity = sqrt(current_data(:,8).^2 + current_data(:,9).^2 + current_data(:,10).^2);
% 
%             % Threshold = 10% of max velocity
%             threshold = 0.10 * max(total_velocity);
% 
%             % Find indices where velocity is above threshold
%             above_threshold = find(total_velocity >= threshold);
% 
%             if isempty(above_threshold)
%                 warning('No samples above threshold for mode %d, trial %d, movement %d', j, i, z);
%                 data_clipped{j,i,z} = current_data; % fallback: keep original
%                 continue
%             end
% 
%             % Clip from first to last index above threshold
%             start_idx = above_threshold(1);
%             end_idx = above_threshold(end);
% 
%             data_clipped{j,i,z} = current_data(start_idx:end_idx, :);
% 
%         end
%     end
% end

%% Clip data based on total velocity threshold
data_clipped = cell(modes, trials, movements);

for j = 1:modes
    for i = 1:trials
        for z = 1:movements
            current_data = data{j,i,z};
            fprintf('%d %d %d\n', j, i,z)
            
            % Compute total velocity magnitude
            total_velocity = sqrt(current_data(:,8).^2 + current_data(:,9).^2 + current_data(:,10).^2);
            
            % Threshold = 10% of max velocity
            threshold = 0.15 * max(total_velocity);
            
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
            
            data_clipped{j,i,z} = current_data(start_idx:end_idx, :);
        end
    end
end


j=1; i=1; z=1;
figure
subplot(2,1,1)
plot(sqrt(data{2,2,2}(:,8).^2+data{2,2,2}(:,9).^2+data{2,2,2}(:,10).^2))
title('Original')
subplot(2,1,2)
plot(sqrt(data_clipped{2,2,2}(:,8).^2+data_clipped{2,2,2}(:,9).^2+data_clipped{2,2,2}(:,10).^2))
title('Clipped')
data = data_clipped;

%% Human effort
groups = [1,2,3,4,5];
grasp_force_mean = zeros(modes,movements);
grasp_force_max = zeros(modes,movements);
grasp_force_std  = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        data_all_trials = [data{j,i,1}(1:end-1,36);data{j,i,2}(1:end-1,36);data{j,i,3}(1:end-1,36);data{j,i,4}(1:end-1,36);data{j,i,5}(1:end-1,36)];
        grasp_force_mean(j,i) = mean(data_all_trials);
        grasp_force_max(j,i) = max(data_all_trials);
        grasp_force_std(j,i)  = std(data_all_trials); 
    end
end


figure()
title('Mean grasp force')
hold on
errorbar(groups, grasp_force_mean(1,:), grasp_force_std(1,:), 'bo-')
errorbar(groups, grasp_force_mean(2,:), grasp_force_std(2,:), 'ro-')
errorbar(groups, grasp_force_mean(3,:), grasp_force_std(3,:), 'go-')
errorbar(groups, grasp_force_mean(4,:), grasp_force_std(4,:), 'ro-.')
errorbar(groups, grasp_force_mean(5,:), grasp_force_std(5,:), 'go-.')
errorbar(groups, grasp_force_mean(6,:), grasp_force_std(6,:), 'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Max grasp force')
hold on
plot(groups, grasp_force_max(1,:),'bo-')
plot(groups, grasp_force_max(2,:),'ro-')
plot(groups, grasp_force_max(3,:),'go-')
plot(groups, grasp_force_max(4,:),'ro-.')
plot(groups, grasp_force_max(5,:),'go-.')
plot(groups, grasp_force_max(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

%% Interaction force
%groups = [1,2,3,4,5,6];
movements=5;

Interaction_force_mean = zeros(modes,movements);
Interaction_force_max = zeros(modes,movements);

Interaction_force_meanX = zeros(modes,movements);
Interaction_force_meanY = zeros(modes,movements);
Interaction_force_meanZ = zeros(modes,movements);

for j = 1:modes
    for i = 1:movements
        data_all_trials = [sqrt(data{j,i,1}(1:end-1,14).^2+data{j,i,1}(1:end-1,15).^2+data{j,i,1}(1:end-1,16).^2); ...
            sqrt(data{j,i,2}(1:end-1,14).^2+data{j,i,2}(1:end-1,15).^2+data{j,i,2}(1:end-1,16).^2); ...
            sqrt(data{j,i,3}(1:end-1,14).^2+data{j,i,3}(1:end-1,15).^2+data{j,i,3}(1:end-1,16).^2); ...
            sqrt(data{j,i,4}(1:end-1,14).^2+data{j,i,4}(1:end-1,15).^2+data{j,i,4}(1:end-1,16).^2); ...
            sqrt(data{j,i,5}(1:end-1,14).^2+data{j,i,5}(1:end-1,15).^2+data{j,i,5}(1:end-1,16).^2)];
        Interaction_force_mean(j,i) = mean(data_all_trials);
        Interaction_force_max(j,i) = max(data_all_trials);

        Interaction_force_meanX(j,i) = mean([abs(data{j,i,1}(1:end-1,14));abs(data{j,i,2}(1:end-1,14)); ...
            abs(data{j,i,3}(1:end-1,14));abs(data{j,i,4}(1:end-1,14));abs(data{j,i,5}(1:end-1,14))]);
        Interaction_force_meanY(j,i) = mean([abs(data{j,i,1}(1:end-1,15));abs(data{j,i,2}(1:end-1,15)); ...
            abs(data{j,i,3}(1:end-1,15));abs(data{j,i,4}(1:end-1,15));abs(data{j,i,5}(1:end-1,15))]);
        Interaction_force_meanZ(j,i) = mean([abs(data{j,i,1}(1:end-1,16));abs(data{j,i,2}(1:end-1,16)); ...
            abs(data{j,i,3}(1:end-1,16));abs(data{j,i,4}(1:end-1,16));abs(data{j,i,5}(1:end-1,16))]);
        % 
        % Interaction_force_meanX(j,i) = mean([abs(data{j,i,1}(1:end-1,14));abs(data{j,i,2}(1:end-1,14)); ...
        %     abs(data{j,i,3}(1:end-1,14));abs(data{j,i,4}(1:end-1,14));abs(data{j,i,5}(1:end-1,14))]);
        % Interaction_force_meanY(j,i) = mean([abs(data{j,i,1}(1:end-1,15));abs(data{j,i,2}(1:end-1,15)); ...
        %     abs(data{j,i,3}(1:end-1,15));abs(data{j,i,4}(1:end-1,15));abs(data{j,i,5}(1:end-1,15))]);
        % Interaction_force_meanZ(j,i) = mean([abs(data{j,i,1}(1:end-1,16));abs(data{j,i,2}(1:end-1,16)); ...
        %     abs(data{j,i,3}(1:end-1,16));abs(data{j,i,4}(1:end-1,16));abs(data{j,i,5}(1:end-1,16))]);
        
    end
end


figure()
title('Mean interaction force')
hold on
plot(groups, Interaction_force_mean(1,:),'bo-')
plot(groups, Interaction_force_mean(2,:),'ro-')
plot(groups, Interaction_force_mean(3,:),'go-')
plot(groups, Interaction_force_mean(4,:),'ro-.')
plot(groups, Interaction_force_mean(5,:),'go-.')
plot(groups, Interaction_force_mean(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Max interaction force')
hold on
plot(groups, Interaction_force_max(1,:),'bo-')
plot(groups, Interaction_force_max(2,:),'ro-')
plot(groups, Interaction_force_max(3,:),'go-')
plot(groups, Interaction_force_max(4,:),'ro-.')
plot(groups, Interaction_force_max(5,:),'go-.')
plot(groups, Interaction_force_max(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Mean interaction force X')
hold on
plot(groups, Interaction_force_meanX(1,:),'bo-')
plot(groups, Interaction_force_meanX(2,:),'ro-')
plot(groups, Interaction_force_meanX(3,:),'go-')
plot(groups, Interaction_force_meanX(4,:),'ro-.')
plot(groups, Interaction_force_meanX(5,:),'go-.')
plot(groups, Interaction_force_meanX(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Mean interaction force Y')
hold on
plot(groups, Interaction_force_meanY(1,:),'bo-')
plot(groups, Interaction_force_meanY(2,:),'ro-')
plot(groups, Interaction_force_meanY(3,:),'go-')
plot(groups, Interaction_force_meanY(4,:),'ro-.')
plot(groups, Interaction_force_meanY(5,:),'go-.')
plot(groups, Interaction_force_meanY(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Mean interaction force Z')
hold on
plot(groups, Interaction_force_meanZ(1,:),'bo-')
plot(groups, Interaction_force_meanZ(2,:),'ro-')
plot(groups, Interaction_force_meanZ(3,:),'go-')
plot(groups, Interaction_force_meanZ(4,:),'ro-.')
plot(groups, Interaction_force_meanZ(5,:),'go-.')
plot(groups, Interaction_force_meanZ(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")



%% Time length
%groups = [1,2,3,4,5,6];
reaching_time = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        data_all_trials = [data{j,i,1}(end,1);data{j,i,2}(end,1);data{j,i,3}(end,1);data{j,i,4}(end,1);data{j,i,5}(end,1)];
        reaching_time(j,i) = mean(data_all_trials);
    end
end



figure()
title('Time to reach')
hold on
plot(groups, reaching_time(1,:),'bo-')
plot(groups, reaching_time(2,:),'ro-')
plot(groups, reaching_time(3,:),'go-')
plot(groups, reaching_time(4,:),'ro-.')
plot(groups, reaching_time(5,:),'go-.')
plot(groups, reaching_time(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

%%
final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(grasp_force_mean(i,:));
    final_results(i,2) = mean(grasp_force_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2"],final_results)
legend("Grasp force across targets mean", "Grasp force across targets max")

final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(Interaction_force_mean(i,:));
    final_results(i,2) = mean(Interaction_force_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2"],final_results)
legend("Interaction force across targets mean", "Interaction force across targets max")

final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(Interaction_force_meanX(i,:));
    final_results(i,2) = mean(Interaction_force_meanY(i,:));
    final_results(i,3) = mean(Interaction_force_meanZ(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2"],final_results)
legend("Interaction force X across targets mean", "Interaction force Y across targets mean",...
    "Interaction force Z across targets mean")

final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(reaching_time(i,:));
    %final_results(i,2) = mean(Interaction_force_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2"],final_results)
legend("Reaching time across targets mean")

%%
figure()
hold on
for i = 1:modes
    plot([i,i,i,i,i], reaching_time(i,:),'o')
end
title('Reaching times across targets')
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")
hold off

figure()
hold on
for i = 1:modes
    plot([i,i,i,i,i], grasp_force_mean(i,:),'o')
end
title('Grasping force across targets')
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")
hold off

%% User data
perceived_effort_comfort = [5,2,5,4,2,3; 6,7,7,7,7,7]; % Subject 21
%perceived_effort_comfort = [6,7,7,5,4; 4,4,3,4,5]; % Subject 9
%perceived_effort_comfort = [3,1,2,3,2; 6,7,6,5,7]; % Subject 8
%comfort = [4,4,5,7,7];
figure()
title('Perceived effort')
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed", "Admittance 2"],transpose(perceived_effort_comfort))
legend("Perceived effort", "Comfort")

% figure()
% title('Perceived effort')
% bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],transpose(comfort))
% legend("Comfort")


%% Current
%groups = [1,2,3,4,5,6];
Current_mean = zeros(modes,movements);
Current_max = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        d1 = [0];
        d2 = [0];
        d3 = [0];
        d4 = [0];
        d5 = [0];
        for z = 37:43
            d1 = d1 +data{j,i,1}(1:end-1,z).^2;
        d2 = d2 +data{j,i,2}(1:end-1,z).^2;
        d3 = d3 +data{j,i,3}(1:end-1,z).^2;
        d4 = d4 + data{j,i,4}(1:end-1,z).^2;
        d5 = d5 + data{j,i,5}(1:end-1,z).^2;
        end

        data_all_trials = [sqrt(d1);sqrt(d2);sqrt(d3);sqrt(d4);sqrt(d5)];
        Current_mean(j,i) = mean(data_all_trials);
        Current_max(j,i) = max(data_all_trials);
        
    end
end


figure()
title('Mean current')
hold on
plot(groups, Current_mean(1,:),'bo-')
plot(groups, Current_mean(2,:),'ro-')
plot(groups, Current_mean(3,:),'go-')
plot(groups, Current_mean(4,:),'ro-.')
plot(groups, Current_mean(5,:),'go-.')
plot(groups, Current_mean(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")

figure()
title('Max current')
hold on
plot(groups, Current_max(1,:),'bo-')
plot(groups, Current_max(2,:),'ro-')
plot(groups, Current_max(3,:),'go-')
plot(groups, Current_max(4,:),'ro-.')
plot(groups, Current_max(5,:),'go-.')
plot(groups, Current_max(6,:),'bo-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2")


final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(Current_mean(i,:));
    final_results(i,2) = mean(Current_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed","Admittance 2"],final_results)
legend("Current across targets mean", "Current across targets max")


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
% 
% %% Calculate the error in prediction position
% time = data{i,1}(:,1);
% predErrorX = zeros(length(time),1);
% predErrorY = zeros(length(time),1);
% predErrorZ = zeros(length(time),1);
% for j = 1:length(time)
%     predX = data{i,1}(j,20);
%     predY = data{i,1}(j,21);
%     predZ = data{i,1}(j,22);
%     s_diff = data{i,1}(j,32);
% 
%     [~, idx] = min(abs(time - (time(j)+s_diff)));
%     truePosX = data{i,1}(idx,2);
%     truePosY = data{i,1}(idx,3);
%     truePosZ = data{i,1}(idx,4);
%     predErrorX(j) = abs(predX-truePosX);
%     predErrorY(j) = abs(predY-truePosY); 
%     predErrorZ(j) = abs(predZ-truePosZ); 
% end
% 
% figure()
% subplot(3,1,1)
% plot(time(2:end),predErrorX(2:end)*100)
% title('Pred error for X in cm')
% 
% subplot(3,1,2)
% plot(time(2:end),predErrorY(2:end)*100)
% title('Pred error for Y in cm')
% 
% subplot(3,1,3)
% plot(time(2:end),predErrorZ(2:end)*100)
% title('Pred error for Z in cm')
% 
% %% Calculate the error in prediction Velocity
% time = data{i,1}(:,1);
% predErrorX = zeros(length(time),1);
% predErrorY = zeros(length(time),1);
% predErrorZ = zeros(length(time),1);
% for j = 1:length(time)
%     predX = data{i,1}(j,26);
%     predY = data{i,1}(j,27);
%     predZ = data{i,1}(j,28);
%     s_diff = data{i,1}(j,32);
% 
%     [~, idx] = min(abs(time - (time(j)+s_diff)));
%     trueVelX = data{i,1}(idx,8);
%     trueVelY = data{i,1}(idx,9);
%     trueVelZ = data{i,1}(idx,10);
%     predErrorX(j) = abs(predX-trueVelX);
%     predErrorY(j) = abs(predY-trueVelY); 
%     predErrorZ(j) = abs(predZ-trueVelZ); 
% end
% 
% figure()
% subplot(3,1,1)
% plot(time(2:end),predErrorX(2:end)*100)
% title('Pred velocity error for X in cm')
% 
% subplot(3,1,2)
% plot(time(2:end),predErrorY(2:end)*100)
% title('Pred velocity for Y in cm')
% 
% subplot(3,1,3)
% plot(time(2:end),predErrorZ(2:end)*100)
% title('Pred velocity for Z in cm')
% 
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