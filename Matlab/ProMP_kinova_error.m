
A = load("ProMp/result17.csv");

labels = unique(A(:,1));      % Get unique labels
groups = cell(length(labels),1);
baseData = load("DataCollection/trial2/p1.csv");
%figure()

%%
fig = figure('Name', 'ProMP');
vid = VideoWriter('ProMP.mp4', 'MPEG-4');
vid.FrameRate = 30;
open(vid);

%hold on
horizon = zeros(3,length(labels));
error = zeros(3,length(labels));
error_int = zeros(3,length(labels));
for k = 1:length(labels)
    groups{k} = A(A(:,1) == labels(k), :);
    baseDataCut = baseData(baseData(:,1)>=groups{k}(1,2), :);
    proMP_interp = spline(groups{k}(:,2),groups{k}(:,3),baseDataCut(:,1));
    clf;
    plot(groups{k}(:,2),groups{k}(:,3),'b',baseData(:,1),baseData(:,2),'r',...
        groups{k}(:,2),groups{k}(:,3)+groups{k}(:,4),'b.',...
        groups{k}(:,2),groups{k}(:,3)-groups{k}(:,4),'b.')
    %plot(baseData(:,1),baseData(:,2),'r')
    ylim([0.65, 0.85]);
    xlim([0.0, 3.2]);

    i=1;
    if(length(proMP_interp(:,1))>1)
    while  abs(proMP_interp(i,1)-baseDataCut(i,2))<0.1 %10cm
        horizon(1,k)=baseDataCut(i,1)-baseDataCut(1,1);
        if abs(proMP_interp(i,1)-baseDataCut(i,2))<0.05 %5cm
            horizon(2,k)=baseDataCut(i,1)-baseDataCut(1,1);
            if abs(proMP_interp(i,1)-baseDataCut(i,2))<0.01 %1cm
                horizon(3,k)=baseDataCut(i,1)-baseDataCut(1,1);
            end
        end
        
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end

    i=1;
    if(length(proMP_interp(:,1))>1)
    while  (baseDataCut(i,1)-baseDataCut(1,1))<=0.5 %500ms
        error(1,k)=abs(proMP_interp(i,1)-baseDataCut(i,2));
        error_int(1,k) = error_int(1,k) + abs(proMP_interp(i,1)-baseDataCut(i,2));
        if (baseDataCut(i,1)-baseDataCut(1,1))<=0.1 %100ms
            error(2,k)=abs(proMP_interp(i,1)-baseDataCut(i,2));
            error_int(2,k) = error_int(2,k) + abs(proMP_interp(i,1)-baseDataCut(i,2));
            if (baseDataCut(i,1)-baseDataCut(1,1))<=0.05 %50ms
                error(3,k)=abs(proMP_interp(i,1)-baseDataCut(i,2));
                error_int(3,k) = error_int(3,k) + abs(proMP_interp(i,1)-baseDataCut(i,2));
            end
        end
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end

    drawnow;
    writeVideo(vid,getframe(fig));
end
%hold off
close(vid);

figure()
plot(baseData(1:length(labels),1),horizon(1,:),'b',baseData(1:length(labels),1),baseData(end,1)-baseData(1:length(labels),1),'k',...
    baseData(1:length(labels),1),horizon(2,:),'g',...)
    baseData(1:length(labels),1),horizon(3,:),'r')
legend('Horizon 10cm','Maximum possible horizon','Horizon 5cm','Horizon 1cm')
ylabel('Time horizon (s)')
xlabel('Time of the trajectory (s)')
title('X Horizon until the prediction diverges from the ground truth')

figure()
plot(baseData(1:length(labels),1),error(1,:),'b',...
    baseData(1:length(labels),1),error(2,:),'g',...)
    baseData(1:length(labels),1),error(3,:),'r')
legend('Error at 500ms','Error at 100ms','Error at 50ms')
ylabel('Error (m)')
xlabel('Time of the trajectory (s)')
title('X Error in the future prediction at each timestep')

figure()
plot(baseData(1:length(labels),1),error_int(1,:),'b',...
    baseData(1:length(labels),1),error_int(2,:),'g',...)
    baseData(1:length(labels),1),error_int(3,:),'r')
legend('Integrated error at 500ms','Integrated error at 100ms','Integrated error at 50ms')
ylabel('Integrated error (m)')
xlabel('Time of the trajectory (s)')
title('X Integrated error in the future prediction at each timestep')

disp('X axis')
disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Horizon 10 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(1,1:1063)), max(horizon(1,:)), mean(horizon(1,:)));
fprintf('Horizon 5 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(2,1:1063)), max(horizon(2,:)), mean(horizon(2,:)));
fprintf('Horizon 1 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(3,1:1063)), max(horizon(3,:)), mean(horizon(3,:)));

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error 500ms: min %.4f / max %.4f / mean %.4f \n', min(error(1,1:1063))*100, max(error(1,:))*100, mean(error(1,:))*100);
fprintf('Error 100ms: min %.4f / max %.4f / mean %.4f \n', min(error(2,1:1063))*100, max(error(2,:))*100, mean(error(2,:))*100);
fprintf('Error 50ms: min %.4f / max %.4f / mean %.4f \n', min(error(3,1:1063))*100, max(error(3,:))*100, mean(error(3,:))*100);

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error_int 500ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(1,1:1063))*100, max(error_int(1,:))*100, mean(error_int(1,:))*100);
fprintf('Error_int 100ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(2,1:1063))*100, max(error_int(2,:))*100, mean(error_int(2,:))*100);
fprintf('Error_int 50ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(3,1:1063))*100, max(error_int(3,:))*100, mean(error_int(3,:))*100);

%%
fig = figure('Name', 'ProMP2');
vid = VideoWriter('ProMP2.mp4', 'MPEG-4');
vid.FrameRate = 30;
open(vid);
%baseData = load("DataCollection/trial2/p1.csv");
%figure()
%hold on
horizon = zeros(3,length(labels));
error = zeros(3,length(labels));
error_int = zeros(3,length(labels));
for k = 1:(length(labels)-1)
    groups{k} = A(A(:,1) == labels(k), :);
    baseDataCut = baseData(baseData(:,1)>=groups{k}(1,2), :);
    proMP_interp = spline(groups{k}(:,2),groups{k}(:,7),baseDataCut(:,1));
    clf;
    plot(groups{k}(:,2),groups{k}(:,7),'b',baseData(:,1),baseData(:,3),'r')
    xlim([0.0, 3.2]);

    i=1;
    if(length(proMP_interp(:,1))>1)
    while  abs(proMP_interp(i,1)-baseDataCut(i,3))<0.1 %10cm
        horizon(1,k)=baseDataCut(i,1)-baseDataCut(1,1);
        if abs(proMP_interp(i,1)-baseDataCut(i,3))<0.05 %5cm
            horizon(2,k)=baseDataCut(i,1)-baseDataCut(1,1);
            if abs(proMP_interp(i,1)-baseDataCut(i,3))<0.01 %1cm
                horizon(3,k)=baseDataCut(i,1)-baseDataCut(1,1);
            end
        end
        
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end

    % i=1;
    % if(length(proMP_interp(:,1))>1)
    % while  abs(proMP_interp(i,1)-baseDataCut(i,3))<0.05 %5cm
    %     horizon(2,k)=baseDataCut(i,1)-baseDataCut(1,1);
    %     i=i+1;
    %     if(i > length(proMP_interp(:,1)))
    %         break
    %     end
    % end
    % end
    % 
    % i=1;
    % if(length(proMP_interp(:,1))>1)
    % while  abs(proMP_interp(i,1)-baseDataCut(i,3))<0.01 %1cm
    %     horizon(3,k)=baseDataCut(i,1)-baseDataCut(1,1);
    %     i=i+1;
    %     if(i > length(proMP_interp(:,1)))
    %         break
    %     end
    % end
    % end


    i=1;
    if(length(proMP_interp(:,1))>1)
    while  (baseDataCut(i,1)-baseDataCut(1,1))<=0.5 %500ms
        error(1,k)=abs(proMP_interp(i,1)-baseDataCut(i,3));
        error_int(1,k) = error_int(1,k) + abs(proMP_interp(i,1)-baseDataCut(i,3));
        if (baseDataCut(i,1)-baseDataCut(1,1))<=0.1 %100ms
            error(2,k)=abs(proMP_interp(i,1)-baseDataCut(i,3));
            error_int(2,k) = error_int(2,k) + abs(proMP_interp(i,1)-baseDataCut(i,3));
            if (baseDataCut(i,1)-baseDataCut(1,1))<=0.05 %50ms
                error(3,k)=abs(proMP_interp(i,1)-baseDataCut(i,3));
                error_int(3,k) = error_int(3,k) + abs(proMP_interp(i,1)-baseDataCut(i,3));
            end
        end
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end

    drawnow;
    writeVideo(vid,getframe(fig));
end
%hold off
close(vid);
figure()
plot(baseData(1:length(labels),1),horizon(1,:),'b',baseData(1:length(labels),1),baseData(end,1)-baseData(1:length(labels),1),'k',...
    baseData(1:length(labels),1),horizon(2,:),'g',...)
    baseData(1:length(labels),1),horizon(3,:),'r')
legend('Horizon 10cm','Maximum possible horizon','Horizon 5cm','Horizon 1cm')
ylabel('Time horizon (s)')
xlabel('Time of the trajectory (s)')
title('Y Horizon until the prediction diverges from the ground truth')

figure()
plot(baseData(1:length(labels),1),error(1,:),'b',...
    baseData(1:length(labels),1),error(2,:),'g',...)
    baseData(1:length(labels),1),error(3,:),'r')
legend('Error at 500ms','Error at 100ms','Error at 50ms')
ylabel('Error (m)')
xlabel('Time of the trajectory (s)')
title('Y Error in the future prediction at each timestep')

figure()
plot(baseData(1:length(labels),1),error_int(1,:),'b',...
    baseData(1:length(labels),1),error_int(2,:),'g',...)
    baseData(1:length(labels),1),error_int(3,:),'r')
legend('Integrated error at 500ms','Integrated error at 100ms','Integrated error at 50ms')
ylabel('Integrated error (m)')
xlabel('Time of the trajectory (s)')
title('Y Integrated error in the future prediction at each timestep')

disp('Y axis')
disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Horizon 10 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(1,1:1063)), max(horizon(1,:)), mean(horizon(1,:)));
fprintf('Horizon 5 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(2,1:1063)), max(horizon(2,:)), mean(horizon(2,:)));
fprintf('Horizon 1 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(3,1:1063)), max(horizon(3,:)), mean(horizon(3,:)));

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error 500ms: min %.4f / max %.4f / mean %.4f \n', min(error(1,1:1063))*100, max(error(1,:))*100, mean(error(1,:))*100);
fprintf('Error 100ms: min %.4f / max %.4f / mean %.4f \n', min(error(2,1:1063))*100, max(error(2,:))*100, mean(error(2,:))*100);
fprintf('Error 50ms: min %.4f / max %.4f / mean %.4f \n', min(error(3,1:1063))*100, max(error(3,:))*100, mean(error(3,:))*100);

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error_int 500ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(1,1:1063))*100, max(error_int(1,:))*100, mean(error_int(1,:))*100);
fprintf('Error_int 100ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(2,1:1063))*100, max(error_int(2,:))*100, mean(error_int(2,:))*100);
fprintf('Error_int 50ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(3,1:1063))*100, max(error_int(3,:))*100, mean(error_int(3,:))*100);

%%
fig = figure('Name', 'ProMP3');
vid = VideoWriter('ProMP3.mp4', 'MPEG-4');
vid.FrameRate = 30;
open(vid);
%baseData = load("DataCollection/trial2/p1.csv");
%figure()
%hold on
horizon = zeros(3,length(labels));
error = zeros(3,length(labels));
error_int = zeros(3,length(labels));
for k = 1:length(labels)
    groups{k} = A(A(:,1) == labels(k), :);
    baseDataCut = baseData(baseData(:,1)>=groups{k}(1,2), :);
    proMP_interp = spline(groups{k}(:,2),groups{k}(:,11),baseDataCut(:,1));
    clf;
    plot(groups{k}(:,2),groups{k}(:,11),'b',baseData(:,1),baseData(:,4),'r')
    %plot(baseData(:,1),baseData(:,2),'r')
    %ylim([0.65, 0.85]);
    xlim([0.0, 3.2]);

    i=1;
    if(length(proMP_interp(:,1))>1)
    while  abs(proMP_interp(i,1)-baseDataCut(i,4))<0.1 %10cm
        horizon(1,k)=baseDataCut(i,1)-baseDataCut(1,1);
        if abs(proMP_interp(i,1)-baseDataCut(i,4))<0.05 %5cm
            horizon(2,k)=baseDataCut(i,1)-baseDataCut(1,1);
            if abs(proMP_interp(i,1)-baseDataCut(i,4))<0.01 %1cm
                horizon(3,k)=baseDataCut(i,1)-baseDataCut(1,1);
            end
        end
        
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end

    i=1;
    if(length(proMP_interp(:,1))>1)
    while  (baseDataCut(i,1)-baseDataCut(1,1))<=0.5 %500ms
        error(1,k)=abs(proMP_interp(i,1)-baseDataCut(i,4));
        error_int(1,k) = error_int(1,k) + abs(proMP_interp(i,1)-baseDataCut(i,4));
        if (baseDataCut(i,1)-baseDataCut(1,1))<=0.1 %100ms
            error(2,k)=abs(proMP_interp(i,1)-baseDataCut(i,4));
            error_int(2,k) = error_int(2,k) + abs(proMP_interp(i,1)-baseDataCut(i,4));
            if (baseDataCut(i,1)-baseDataCut(1,1))<=0.05 %50ms
                error(3,k)=abs(proMP_interp(i,1)-baseDataCut(i,4));
                error_int(3,k) = error_int(3,k) + abs(proMP_interp(i,1)-baseDataCut(i,4));
            end
        end
        i=i+1;
        if(i > length(proMP_interp(:,1)))
            break
        end
    end
    end


    drawnow;
    writeVideo(vid,getframe(fig));
end
%hold off
close(vid);

figure()
plot(baseData(1:length(labels),1),horizon(1,:),'b',baseData(1:length(labels),1),baseData(end,1)-baseData(1:length(labels),1),'k',...
    baseData(1:length(labels),1),horizon(2,:),'g',...)
    baseData(1:length(labels),1),horizon(3,:),'r')
legend('Horizon 10cm','Maximum possible horizon','Horizon 5cm','Horizon 1cm')
ylabel('Time horizon (s)')
xlabel('Time of the trajectory (s)')
title('Z Horizon until the prediction diverges from the ground truth')

figure()
plot(baseData(1:length(labels),1),error(1,:),'b',...
    baseData(1:length(labels),1),error(2,:),'g',...)
    baseData(1:length(labels),1),error(3,:),'r')
legend('Error at 500ms','Error at 100ms','Error at 50ms')
ylabel('Error (m)')
xlabel('Time of the trajectory (s)')
title('Z Error in the future prediction at each timestep')

figure()
plot(baseData(1:length(labels),1),error_int(1,:),'b',...
    baseData(1:length(labels),1),error_int(2,:),'g',...)
    baseData(1:length(labels),1),error_int(3,:),'r')
legend('Integrated error at 500ms','Integrated error at 100ms','Integrated error at 50ms')
ylabel('Integrated error (m)')
xlabel('Time of the trajectory (s)')
title('Z Integrated error in the future prediction at each timestep')

disp('Z axis')
disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Horizon 10 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(1,1:1063)), max(horizon(1,:)), mean(horizon(1,:)));
fprintf('Horizon 5 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(2,1:1063)), max(horizon(2,:)), mean(horizon(2,:)));
fprintf('Horizon 1 cm: min %.2f / max %.2f / mean %.2f \n', min(horizon(3,1:1063)), max(horizon(3,:)), mean(horizon(3,:)));

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error 500ms: min %.4f / max %.4f / mean %.4f \n', min(error(1,1:1063))*100, max(error(1,:))*100, mean(error(1,:))*100);
fprintf('Error 100ms: min %.4f / max %.4f / mean %.4f \n', min(error(2,1:1063))*100, max(error(2,:))*100, mean(error(2,:))*100);
fprintf('Error 50ms: min %.4f / max %.4f / mean %.4f \n', min(error(3,1:1063))*100, max(error(3,:))*100, mean(error(3,:))*100);

disp('We discard for the min the last 0.5 second, since if not it would always be 0')
fprintf('Error_int 500ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(1,1:1063))*100, max(error_int(1,:))*100, mean(error_int(1,:))*100);
fprintf('Error_int 100ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(2,1:1063))*100, max(error_int(2,:))*100, mean(error_int(2,:))*100);
fprintf('Error_int 50ms: min %.4f / max %.4f / mean %.4f \n', min(error_int(3,1:1063))*100, max(error_int(3,:))*100, mean(error_int(3,:))*100);