% data {mode, movement}
modes = 5;
movements = 1;

data = cell(5,1);
for j = 1:modes
    for i = 1:movements
        data{j,i} = load(sprintf("tests/subject10mode%dmov%d.csv", j,i));
        %fprintf('%d %d\n', j, i)
    end
end


% %% admittance data
% admittance_data = cell(4,5);
% for i = 1:movements
%     admittance_data{1,i} = data{1,i};
%     %fprintf('%d %d\n', j, i)
% end
% for j = [3,4,5]
%     for i = 1:movements
%         admittance_data{j-1,i} = readmatrix(sprintf("subject5/subject5mode1mov%d%d.csv", j,i));
%         %fprintf('%d %d\n', j, i)
%     end
% end

%% Human effort
groups = [1,2,3,4,5];
grasp_force_mean = zeros(modes,movements);
grasp_force_max = zeros(modes,movements);
grasp_force_std  = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        grasp_force_mean(j,i) = mean(data{j,i}(:,36));
        grasp_force_max(j,i) = max(data{j,i}(:,36));
        grasp_force_std(j,i)  = std(data{j,i}(:,36)); 
    end
end

% for i = 1:movements
%     grasp_force_mean(1,i) = mean([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_max(1,i) = max([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_std(1,i)  = std([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
% end


figure()
title('Mean grasp force')
hold on
errorbar(groups, grasp_force_mean(1,:), grasp_force_std(1,:), 'bo-')
errorbar(groups, grasp_force_mean(2,:), grasp_force_std(2,:), 'ro-')
errorbar(groups, grasp_force_mean(3,:), grasp_force_std(3,:), 'go-')
errorbar(groups, grasp_force_mean(4,:), grasp_force_std(4,:), 'ro-.')
errorbar(groups, grasp_force_mean(5,:), grasp_force_std(5,:), 'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

figure()
title('Max grasp force')
hold on
plot(groups, grasp_force_max(1,:),'bo-')
plot(groups, grasp_force_max(2,:),'ro-')
plot(groups, grasp_force_max(3,:),'go-')
plot(groups, grasp_force_max(4,:),'ro-.')
plot(groups, grasp_force_max(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

%% Interaction force
groups = [1,2,3,4,5];
Interaction_force_mean = zeros(modes,movements);
Interaction_force_max = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        Interaction_force_mean(j,i) = mean(sqrt(data{j,i}(:,14).^2+data{j,i}(:,15).^2+data{j,i}(:,16).^2));
        Interaction_force_max(j,i) = max(sqrt(data{j,i}(:,14).^2+data{j,i}(:,15).^2+data{j,i}(:,16).^2));
        
    end
end

% for i = 1:movements
%     Interaction_force_mean(1,i) = mean([sqrt(admittance_data{1,i}(:,14).^2+admittance_data{1,i}(:,15).^2+admittance_data{1,i}(:,16).^2); ...
%         sqrt(admittance_data{2,i}(:,14).^2+admittance_data{2,i}(:,15).^2+admittance_data{2,i}(:,16).^2);...
%         sqrt(admittance_data{3,i}(:,14).^2+admittance_data{3,i}(:,15).^2+admittance_data{3,i}(:,16).^2);...
%         sqrt(admittance_data{4,i}(:,14).^2+admittance_data{4,i}(:,15).^2+admittance_data{4,i}(:,16).^2);]);
%     Interaction_force_max(1,i) = max([sqrt(admittance_data{1,i}(:,14).^2+admittance_data{1,i}(:,15).^2+admittance_data{1,i}(:,16).^2); ...
%         sqrt(admittance_data{2,i}(:,14).^2+admittance_data{2,i}(:,15).^2+admittance_data{2,i}(:,16).^2);...
%         sqrt(admittance_data{3,i}(:,14).^2+admittance_data{3,i}(:,15).^2+admittance_data{3,i}(:,16).^2);...
%         sqrt(admittance_data{4,i}(:,14).^2+admittance_data{4,i}(:,15).^2+admittance_data{4,i}(:,16).^2);]);
% end



figure()
title('Mean interaction force')
hold on
plot(groups, Interaction_force_mean(1,:),'bo-')
plot(groups, Interaction_force_mean(2,:),'ro-')
plot(groups, Interaction_force_mean(3,:),'go-')
plot(groups, Interaction_force_mean(4,:),'ro-.')
plot(groups, Interaction_force_mean(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

figure()
title('Max interaction force')
hold on
plot(groups, Interaction_force_max(1,:),'bo-')
plot(groups, Interaction_force_max(2,:),'ro-')
plot(groups, Interaction_force_max(3,:),'go-')
plot(groups, Interaction_force_max(4,:),'ro-.')
plot(groups, Interaction_force_max(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

%% Time length
groups = [1,2,3,4,5];
reaching_time = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        reaching_time(j,i) = data{j,i}(end,1);
    end
end

% for i = 1:movements
%     reaching_time(1,i) = mean([admittance_data{1,i}(end,1);admittance_data{2,i}(end,1);admittance_data{3,i}(end,1);admittance_data{4,i}(end,1)]);
% end


figure()
title('Time to reach')
hold on
plot(groups, reaching_time(1,:),'bo-')
plot(groups, reaching_time(2,:),'ro-')
plot(groups, reaching_time(3,:),'go-')
plot(groups, reaching_time(4,:),'ro-.')
plot(groups, reaching_time(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

%%
final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(grasp_force_mean(i,:));
    final_results(i,2) = mean(grasp_force_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],final_results)
legend("Grasp force across targets mean", "Grasp force across targets max")

final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(Interaction_force_mean(i,:));
    final_results(i,2) = mean(Interaction_force_max(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],final_results)
legend("Interaction force across targets mean", "Interaction force across targets max")


%%
modes = 5;
movements = 5;

data_obstacle = cell(5,5);
for j = 1:modes
    for i = 1:movements
        data_obstacle{j,i} = load(sprintf("subject5/subject5mode%dmov%d.csv", j,i));
        if (j>1)
            data_obstacle{j,i} = readmatrix(sprintf("subject5/subject5mode%dmov10%d.csv", j,i));
        end
        %fprintf('%d %d\n', j, i)
    end
end

%% Human effort
groups = [1,2,3,4,5];
grasp_force_mean = zeros(modes,movements);
grasp_force_max = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        %fprintf('%d %d\n', j, i)
        grasp_force_mean(j,i) = mean(data_obstacle{j,i}(1:end-1,36));
        grasp_force_max(j,i) = max(data_obstacle{j,i}(1:end-1,36));
        
    end
end

% for i = 1:movements
%     grasp_force_mean(1,i) = mean([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_max(1,i) = max([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_std(1,i)  = std([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
% end

figure()
title('Mean grasp force')
hold on
plot(groups, grasp_force_mean(1,:),'bo-')
plot(groups, grasp_force_mean(2,:),'ro-')
plot(groups, grasp_force_mean(3,:),'go-')
plot(groups, grasp_force_mean(4,:),'ro-.')
plot(groups, grasp_force_mean(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

figure()
title('Max grasp force')
hold on
plot(groups, grasp_force_max(1,:),'bo-')
plot(groups, grasp_force_max(2,:),'ro-')
plot(groups, grasp_force_max(3,:),'go-')
plot(groups, grasp_force_max(4,:),'ro-.')
plot(groups, grasp_force_max(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

%% Interaction force
groups = [1,2,3,4,5];
Interaction_force_mean = zeros(modes,movements);
Interaction_force_max = zeros(modes,movements);

for j = 1:modes
    for i = 1:movements
        Interaction_force_mean(j,i) = mean(sqrt(data_obstacle{j,i}(:,14).^2+data_obstacle{j,i}(:,15).^2+data_obstacle{j,i}(:,16).^2));
        Interaction_force_max(j,i) = max(sqrt(data_obstacle{j,i}(:,14).^2+data_obstacle{j,i}(:,15).^2+data_obstacle{j,i}(:,16).^2));

        
    end
end

% for i = 1:movements
%     Interaction_force_mean(1,i) = mean([sqrt(admittance_data{1,i}(:,14).^2+admittance_data{1,i}(:,15).^2+admittance_data{1,i}(:,16).^2); ...
%         sqrt(admittance_data{2,i}(:,14).^2+admittance_data{2,i}(:,15).^2+admittance_data{2,i}(:,16).^2);...
%         sqrt(admittance_data{3,i}(:,14).^2+admittance_data{3,i}(:,15).^2+admittance_data{3,i}(:,16).^2);...
%         sqrt(admittance_data{4,i}(:,14).^2+admittance_data{4,i}(:,15).^2+admittance_data{4,i}(:,16).^2);]);
%     Interaction_force_max(1,i) = max([sqrt(admittance_data{1,i}(:,14).^2+admittance_data{1,i}(:,15).^2+admittance_data{1,i}(:,16).^2); ...
%         sqrt(admittance_data{2,i}(:,14).^2+admittance_data{2,i}(:,15).^2+admittance_data{2,i}(:,16).^2);...
%         sqrt(admittance_data{3,i}(:,14).^2+admittance_data{3,i}(:,15).^2+admittance_data{3,i}(:,16).^2);...
%         sqrt(admittance_data{4,i}(:,14).^2+admittance_data{4,i}(:,15).^2+admittance_data{4,i}(:,16).^2);]);
% end


figure()
title('Mean interaction force')
hold on
plot(groups, Interaction_force_mean(1,:),'bo-')
plot(groups, Interaction_force_mean(2,:),'ro-')
plot(groups, Interaction_force_mean(3,:),'go-')
plot(groups, Interaction_force_mean(4,:),'ro-.')
plot(groups, Interaction_force_mean(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")

figure()
title('Max interaction force')
hold on
plot(groups, Interaction_force_max(1,:),'bo-')
plot(groups, Interaction_force_max(2,:),'ro-')
plot(groups, Interaction_force_max(3,:),'go-')
plot(groups, Interaction_force_max(4,:),'ro-.')
plot(groups, Interaction_force_max(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")



%% Time length
groups = [1,2,3,4,5];
reaching_time = zeros(modes,movements);
Interaction_force_max = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        reaching_time(j,i) = data_obstacle{j,i}(end,1);
    end
end

% for i = 1:movements
%     reaching_time(1,i) = mean([admittance_data{1,i}(end,1);admittance_data{2,i}(end,1);admittance_data{3,i}(end,1);admittance_data{4,i}(end,1)]);
% end

figure()
title('Time to reach')
hold on
plot(groups, reaching_time(1,:),'bo-')
plot(groups, reaching_time(2,:),'ro-')
plot(groups, reaching_time(3,:),'go-')
plot(groups, reaching_time(4,:),'ro-.')
plot(groups, reaching_time(5,:),'go-.')
hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed")


%% User data
perceived_effort = [6,6,6,5,5; 6,5,5,4,6];
comfort = [4,4,4,5,5 ; 4, 4,5,5,3];
figure()
title('Perceived effort')
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],transpose(perceived_effort))
legend("Perceived effort no obstacle", "Perceived effort obstacle")

figure()
title('Perceived effort')
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],transpose(comfort))
legend("Comfort no obstacle", "Comfort obstacle")

%% COMPARISON WITH AND WITHOUT OBSTACLE
groups = [1,2,3,4,5];
grasp_force_mean = zeros(modes,movements);
grasp_force_max = zeros(modes,movements);
grasp_force_mean_obs = zeros(modes,movements);
grasp_force_max_obs = zeros(modes,movements);
for j = 1:modes
    for i = 1:movements
        grasp_force_mean(j,i) = mean(data{j,i}(:,36));
        grasp_force_max(j,i) = max(data{j,i}(:,36));

        grasp_force_mean_obs(j,i) = mean(data_obstacle{j,i}(1:end-1,36));
        grasp_force_max_obs(j,i) = max(data_obstacle{j,i}(1:end-1,36));
        
    end
end

% for i = 1:movements
%     grasp_force_mean(1,i) = mean([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_max(1,i) = max([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_mean_obs(1,i) = mean([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
%     grasp_force_max_obs(1,i) = max([admittance_data{1,i}(1:end-1,36);admittance_data{2,i}(1:end-1,36);admittance_data{3,i}(1:end-1,36);admittance_data{4,i}(1:end-1,36)]);
% end


figure()
title('Mean grasp force')
hold on
plot(groups, grasp_force_mean(1,:),'bo-')
plot(groups, grasp_force_mean(2,:),'ro-')
plot(groups, grasp_force_mean(3,:),'go-')
plot(groups, grasp_force_mean(4,:),'co-')
plot(groups, grasp_force_mean(5,:),'mo-')
plot(groups, grasp_force_mean_obs(1,:),'bo-.')
plot(groups, grasp_force_mean_obs(2,:),'ro-.')
plot(groups, grasp_force_mean_obs(3,:),'go-.')
plot(groups, grasp_force_mean_obs(4,:),'co-.')
plot(groups, grasp_force_mean_obs(5,:),'mo-.')

hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed", ...
    "Obs-Admittance", "Obs-ProMP var", "Obs-Kalman var", "Obs-ProMP fixed", "Obs-Kalman fixed")


figure()
title('Max grasp force')
hold on
plot(groups, grasp_force_max(1,:),'bo-')
plot(groups, grasp_force_max(2,:),'ro-')
plot(groups, grasp_force_max(3,:),'go-')
plot(groups, grasp_force_max(4,:),'co-')
plot(groups, grasp_force_max(5,:),'mo-')
plot(groups, grasp_force_max_obs(1,:),'bo-.')
plot(groups, grasp_force_max_obs(2,:),'ro-.')
plot(groups, grasp_force_max_obs(3,:),'go-.')
plot(groups, grasp_force_max_obs(4,:),'co-.')
plot(groups, grasp_force_max_obs(5,:),'mo-.')

hold off
legend("Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed", ...
    "Obs-Admittance", "Obs-ProMP var", "Obs-Kalman var", "Obs-ProMP fixed", "Obs-Kalman fixed")


final_results = zeros(5,1);
for i = 1:modes
    final_results(i,1) = mean(grasp_force_mean(i,:));
    final_results(i,2) = mean(grasp_force_max(i,:));
    final_results(i,3) = mean(grasp_force_mean_obs(i,:));
    final_results(i,4) = mean(grasp_force_max_obs(i,:));
end
figure()
bar(["Admittance", "ProMP var", "Kalman var", "ProMP fixed", "Kalman fixed"],final_results)
legend("Grasp force across targets mean", "Grasp force across targets max", ...
    "Obs grasp force across targets mean", "Obs grasp force across targets max")

%%
clear

data{1,1} = load("tests/subject12mode3mov1.csv");
data{2,1} = load("var_horizon1/r2.csv");
data{3,1} = load("var_horizon1/r3.csv");

i = 1;
s_diff = 0.030;
%

figure()
subplot(3,1,1)
plot(data{i,1}(:,1),data{i,1}(:,2))
title('posX')

subplot(3,1,2)
plot(data{i,1}(:,1),data{i,1}(:,3))
title('posY')

subplot(3,1,3)
plot(data{i,1}(:,1),data{i,1}(:,4))
title('posZ')

figure()
subplot(3,1,1)
plot(data{i,1}(:,1),data{i,1}(:,8))
title('VelX')

subplot(3,1,2)
plot(data{i,1}(:,1),data{i,1}(:,9))
title('VelY')

subplot(3,1,3)
plot(data{i,1}(:,1),data{i,1}(:,10))
title('VelZ')

figure
subplot(3,1,1)
hold on
plot(data{i,1}(:,1),data{i,1}(:,8)./max(data{i,1}(:,8)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio VelX')
hold off

subplot(3,1,2)
hold on
plot(data{i,1}(:,1),data{i,1}(:,9)./max(data{i,1}(:,9)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio VelY')
hold off

subplot(3,1,3)
hold on
plot(data{i,1}(:,1),data{i,1}(:,10)./max(abs(data{i,1}(:,10))))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(abs(data{i,1}(:,36))))
title('ratio VelZ')
hold off

% figure()
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,2))
% plot(data{i,1}(:,1),data{i,1}(:,3))
% plot(data{i,1}(:,1),data{i,1}(:,4))
% plot(data{i,1}(:,1),data{i,1}(:,5))
% plot(data{i,1}(:,1),data{i,1}(:,6))
% plot(data{i,1}(:,1),data{i,1}(:,7))
% hold off
% legend("posX","posY","posZ","rot1","rot2","rot3")
% title('True trajectory')
% 
% 
% figure()
% hold on
% plot(data{i,1}(:,1),data{i,1}(:,8))
% plot(data{i,1}(:,1),data{i,1}(:,9))
% plot(data{i,1}(:,1),data{i,1}(:,10))
% plot(data{i,1}(:,1),data{i,1}(:,11))
% plot(data{i,1}(:,1),data{i,1}(:,12))
% plot(data{i,1}(:,1),data{i,1}(:,13))
% hold off
% legend("velX","velY","velZ","velrot1","velrot2","velrot3")
% title('True velocity')

figure()
plot(data{i,1}(:,1), sqrt((data{i,1}(:,8).^2+data{i,1}(:,9).^2+data{i,1}(:,10).^2)))
title('Total velocity')


figure()
subplot(3,1,1)
plot(data{i,1}(:,1),data{i,1}(:,14))
title('TorqueX')

subplot(3,1,2)
plot(data{i,1}(:,1),data{i,1}(:,15))
title('TorqueY')

subplot(3,1,3)
plot(data{i,1}(:,1),data{i,1}(:,16)+30)
title('TorqueZ')

figure
subplot(3,1,1)
hold on
plot(data{i,1}(:,1),data{i,1}(:,14)./max(data{i,1}(:,14)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueX')
hold off

subplot(3,1,2)
hold on
plot(data{i,1}(:,1),data{i,1}(:,15)./max(data{i,1}(:,15)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueY')
hold off


subplot(3,1,3)
hold on
plot(data{i,1}(:,1),data{i,1}(:,16)./max(abs(data{i,1}(:,16))))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(abs(data{i,1}(:,36))))
title('ratio TorqueZ')
hold off


figure()
subplot(3,1,1)

plot(data{i,1}(:,1),data{i,1}(:,14).*data{i,1}(:,8))
title('TorqueX*VelocityX')


subplot(3,1,2)
plot(data{i,1}(:,1),data{i,1}(:,15).*data{i,1}(:,9))
title('TorqueY*VelocityY')

subplot(3,1,3)
plot((data{i,1}(:,1)+30),data{i,1}(:,16).*data{i,1}(:,10))
title('TorqueZ*VelocityZ')

figure
subplot(3,1,1)
hold on
plot(data{i,1}(:,1),data{i,1}(:,14).*data{i,1}(:,8)./max(data{i,1}(:,14).*data{i,1}(:,8)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueX*VelocityX')
hold off

subplot(3,1,2)
hold on
plot(data{i,1}(:,1),data{i,1}(:,15).*data{i,1}(:,9)./max(data{i,1}(:,15).*data{i,1}(:,9)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueY*VelocityY')
hold off


subplot(3,1,3)
hold on
plot(data{i,1}(:,1),data{i,1}(:,16).*data{i,1}(:,10)./max(data{i,1}(:,16).*data{i,1}(:,10)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueZ*VelocityZ')
hold off

%
figure()
subplot(3,1,1)
plot(data{i,1}(:,1),data{i,1}(:,33))
title('Kx')

subplot(3,1,2)
plot(data{i,1}(:,1),data{i,1}(:,34))
title('Ky')

subplot(3,1,3)
plot(data{i,1}(:,1),data{i,1}(:,35))
title('Kz')

% sensor reading
figure()
plot(data{i,1}(:,1),data{i,1}(:,36))
title('Human effort')



idx = abs(data{i,1}(:,8)) < 0.1;
figure
subplot(3,1,1)
hold on
plot(data{i,1}(:,1),data{i,1}(:,14).*(idx)./max(data{i,1}(:,14)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueX*(1/VelocityX)')
hold off

idx = abs(data{i,1}(:,9)) < 0.1;
subplot(3,1,2)
hold on
plot(data{i,1}(:,1),data{i,1}(:,15).*idx./max(data{i,1}(:,15)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueY*1/VelocityY')
hold off

idx = abs(data{i,1}(:,10)) < 0.1;
subplot(3,1,3)
hold on
plot(data{i,1}(:,1),data{i,1}(:,16).*idx./max(data{i,1}(:,16)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
title('ratio TorqueZ*1/VelocityZ')
hold off


figure()
hold on
plot(data{i,1}(:,1),data{i,1}(:,8)./max(data{i,1}(:,8)))
plot(data{i,1}(:,1),data{i,1}(:,14)./max(data{i,1}(:,14)))
plot(data{i,1}(:,1),data{i,1}(:,33)./max(data{i,1}(:,33)))
plot(data{i,1}(:,1),data{i,1}(:,36)./max(data{i,1}(:,36)))
legend('vel', 'force','stiffness','effort')
hold off

figure()
hold on
plot(data{i,1}(:,1),data{i,1}(:,2))
plot(data{i,1}(:,1),data{i,1}(:,20))
plot(data{i,1}(:,1),sign(data{i,1}(:,8)))
plot(data{i,1}(:,1),data{i,1}(:,8))

hold off


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

