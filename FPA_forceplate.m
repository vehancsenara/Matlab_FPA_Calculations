function [foot_tangent,side_of_the_emg] = FPA_forceplate(fcp, EM, timetable)

%EMG ACTIVITY -MEAN - STANDARD DEVIATION
    avgg=mean(EM.r_Gastrocnem_);
    stdg=std(EM.r_Gastrocnem_);
    idk=1:length(EM.r_Gastrocnem_);%length of the data
    xx=(1:length(EM.r_Gastrocnem_))/600;
    c=3; % constant
    %AQUIRING THE LOCAL MAXIMAS AND FILTER ONLY THOSE WHICH ARE ABOVE THE
    %THRESHOLD
    locmax_0=islocalmax(EM.r_Gastrocnem_);
    some=EM.r_Gastrocnem_(locmax_0)>avgg+stdg*c;%this is the threshold
    log_ic=idk(locmax_0);
    filter_ix=log_ic(some);
    % PLOT THE THRESHOLD AND THE EMG SIGNAL
    figure
    hold on
    title("SCI - EMG's local maximas above the threshold", fontsize=38, FontWeight='bold')   
    plot(xx,EM.r_Gastrocnem_, xx(filter_ix), EM.r_Gastrocnem_(filter_ix), 'r*',LineWidth=2)
    grid on
    grid minor
    yline(avgg+c*stdg,LineWidth=2, Color="red")
    ylabel("sEMG [mV]", fontsize=26, FontWeight='bold')
    xlabel ("Time [s]", fontsize=26, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([0 7])
     legend(["Right Gastrocnemius ", "Local maximas","Threshold"],  'Location', 'southeastoutside')
    % GET WHICH SIDE THE EMG ELECTRODES WERE
    side_of_the_emg=EM.Properties.VariableNames{1}(1);
    for k=1:3
    lr_logic=ismember(timetable(:,k)*600, filter_ix);
    yy=find(lr_logic==1);
    if sum(yy)==0
        if side_of_the_emg=='r'
            disp("That's the left side")            
        else
            disp("That's the right side-not")
        end
    else
        if side_of_the_emg=='r'
            disp("That's the right side-tab")
        else
            disp("That's the left side other")
        end
    end
    end
 if side_of_the_emg=='r'
     side_of_the_emg='right';
 else
     side_of_the_emg='left';
 end
foot_tangent=zeros(1,3);
for k=1:3
HCx=fcp.(2+k)(1:round(size(fcp.(2+k),1)*0.05));
TOx=fcp.(2+k)((round(size(fcp.(2+k),1)*0.95)):(round(size(fcp.(2+k),1)*1)));
HCy=fcp.(6+k)(1:(round(size(fcp.(6+k),1)*0.05)));
TOy=fcp.(6+k)((round(size(fcp.(6+k),1)*0.95)):(round(size(fcp.(6+k),1)*1)));
  
% some tangent calculation with size adjustment again
if size(HCx,1)== size(TOx,1)
foot_x= HCx-TOx;
foot_y= HCy-TOy;
foot_tangent(1,k)=mean(atand(foot_x./foot_y));

elseif size(HCx,1)> size(TOx,1)
       HCx=HCx(1:size(TOx, 1));
       HCy=HCy(1:size(TOx, 1));
       foot_x= HCx-TOx;
       foot_y= HCy-TOy;
       
       foot_tangent(1,k)=mean(atand(foot_x./foot_y));
elseif size(HCx,1)<size(TOx,1)
       TOx=TOx(1:size(HCx, 1));
       TOy=TOy(1:size(HCx, 1));
       foot_x= HCx-TOx;
       foot_y= HCy-TOy;
       foot_tangent(1,k)=mean(atand(foot_x./foot_y));
end
end
end