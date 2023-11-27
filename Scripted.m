%'p02_r_1stp1fcp_slow_2.xlsx' mashogy vannak cimkezve az emg-k
clear;
close all;
% READ TABLES IN
%SCIS
% f_name='p02_r_1stp1fcp_slow_1.xlsx';
% fcp=readtable(f_name,'Sheet','fcp');
% emg=readtable(f_name,'Sheet','emg');
% kin=readtable(f_name,'Sheet','kin');


%PREFES
% f_name='p02_r_1stp1fcp_norm_1.xlsx';
% fcp=readtable(f_name,'Sheet','fcp');
% emg=readtable(f_name,'Sheet','emg');
% kin=readtable(f_name,'Sheet','kin');

%POSTFES
% 
f_name='post02_l_1spt1fcp_norm_4.xlsx';
fcp=readtable(f_name,'Sheet','fcp');
emg=readtable(f_name,'Sheet','emg');
kin=readtable(f_name,'Sheet','kin');

%HELT
% f_name='h06_r_1stp1fcp_norm_2.xlsx';
% fcp=readtable(f_name,'Sheet','fcp');
% emg=readtable(f_name,'Sheet','emg');
% kin=readtable(f_name,'Sheet','kin');

%FILTERS FOR CERTAIN COLUMNS OF THE TABLE
hasMatchfz = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'F_z', 'once'));
hasMatchcx = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'CoP_x', 'once'));
hasMatchcy = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'CoP_y', 'once'));
hasMatchk = ~cellfun('isempty', regexp(kin.Properties.VariableNames,'HEE|TOE', 'once'));
hasMatche = ~cellfun('isempty', regexp(emg.Properties.VariableNames,'Gastrocnem|Soleus|Vast_Lat', 'once'));
%ACTUALLY FILTERING
Fz=fcp(:, fcp.Properties.VariableNames(hasMatchfz));
Copx=fcp(:, fcp.Properties.VariableNames(hasMatchcx));
Copy=fcp(:, fcp.Properties.VariableNames(hasMatchcy));
kinfoot=kin(:, kin.Properties.VariableNames(hasMatchk));
EM=emg(:, emg.Properties.VariableNames(hasMatche));
[Fz, Copy, Copx]=order_fz(Fz, Copy, Copx, length(EM.r_Gastrocnem_));

%FPA CALCULATION AND ADJUSTING THE COPs
[FPA_left,FPA_right, right_tangent, left_tangent, left, right, upsampled_kin] = FPA_kin(kinfoot, length(EM.r_Gastrocnem_));

Copy.fp1_CoP_y=Copy.fp1_CoP_y-mean(Copy.fp1_CoP_y(1:500));
Copy.fp2_CoP_y=Copy.fp2_CoP_y-mean(Copy.fp2_CoP_y(1:500));
Copy.fp3_CoP_y=Copy.fp3_CoP_y-mean(Copy.fp3_CoP_y(1:500));

% LOWPASS FILTERING THE DATA
d1=designfilt("lowpassiir",FilterOrder=6, ...
    HalfPowerFrequency=0.05,DesignMethod="butter");
d2=designfilt('bandpassiir', ...
    'FilterOrder',4,'HalfPowerFrequency1',20, ...
    'HalfPowerFrequency2',400,'SampleRate',1200);
d3=designfilt("lowpassiir",FilterOrder=6, ...
    HalfPowerFrequency=0.12,DesignMethod="butter");
%LOOP FOR FILTERING
old_em=EM;
for k=1:3
fz=filtfilt(d1, Fz.(k));
cx=filtfilt(d1, Copx.(k));
cy=filtfilt(d1, Copy.(k));
em=filtfilt(d2,EM.(k));
em1=smoothEMG(em, 55);
em2=filtfilt(d3, em1);
EM.(k)=em2/max(em2);
Fz.(k)=fz;
Copx.(k)=cx;
Copy.(k)=cy;
end

tem=(1:1:length(EM.r_Vast_Lat_))/1200;
figure 
hold on
sgtitle("Post-FES - Raw and processed EMG", fontsize=36, FontWeight='bold')
subplot(3,1,1)
hold on
plot(tem, old_em.r_Vast_Lat_, Color='blue')
plot(tem, EM.r_Vast_Lat_, LineWidth=3,  Color='red')
grid on
grid minor
xlim([1/1200 length(EM.r_Vast_Lat_)/1200])
ylabel("sEMG [V]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
     legend(["Right Vastus Lateralis raw","Right Vastus Lateralis processed"],  'Location', 'southeastoutside')
    
subplot(3,1,2)
hold on
plot(tem, old_em.r_Gastrocnem_, Color='blue')
plot(tem, EM.r_Gastrocnem_,LineWidth=3,  Color='red')
grid on
grid minor
ylabel("sEMG [V]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize',22);
    xlim([1/1200 length(EM.r_Gastrocnem_)/1200])
    legend(["Right Gastrocnemius raw", "Right Gastrocnemius processed"], 'Location', 'southeastoutside')
   
subplot(3,1,3)
hold on
plot(tem, old_em.r_Soleus, Color='blue')
plot(tem, EM.r_Soleus,LineWidth=3, Color='red')
grid on
grid minor
xlim([1/1200 length(EM.r_Soleus)/1200])
ylabel("sEMG [V]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    legend(["Right Soleus raw", "Right Soleus processed"],  'Location', 'southeastoutside')
    

% GETTING THE MEAN AND STANDARD DEVIATION
stdf=zeros(3,1);
stdc=zeros(6,1);
meanf=zeros(3,1);
meanc=zeros(6,1);

%LOOP FOR THE MEAN AND STANDARD DEVIATION
for k=1:3
sf=std(Fz.(k)(1:500));
scx=std(Copx.(k)(1:500));
scy=std(Copy.(k)(1:500));
stdf(k,1)=sf;
stdc(k,1)=scx;
stdc(3+k, 1)=scy;
mf=mean(Fz.(k)(1:500));
mcx=mean(Copx.(k)(1:500));
mcy=mean(Copy.(k)(1:500));
meanf(k,1)=mf;
meanc(k,1)=mcx;
meanc(3+k, 1)=mcy;
end
%CONSTANT
c=adapt_c(Fz, stdf,meanf);
%
%PLOTTING THE RESULTS OF THE THRESHOLD WHICH IS THE MEAN
tlab=(1:1:length(Copx.fp1_CoP_x))/600;
fz_name=["Force plate_1 Fz", "Force plate_2 Fz", "Force plate_3 Fz"];
copy_name=["Force plate_1 CoP_y", "Force plate_2 CoP_y", "Force plate_3 CoP_y"];
copx_name=["Force plate_1 CoP_x", "Force plate_2 CoP_x", "Force plate_3 CoP_x"];
    figure
    hold on
    sgtitle("Post-FES - Forceplate Force component z", fontsize=38, FontWeight='bold')
    for k=1:3
    subplot(3,1,k)
    plot(tlab, Fz.(k), LineWidth=3)
    grid on
grid minor
    ylabel("Force [N]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
   % xlim([1 6]) %helt
    xlim([2 13]) %scis
   % xlim([2 11]) %pre
    %xlim([4 10.5]) %post
    yline(meanf(k,1)-c(k,1)*stdf(k,1),LineWidth=2, Color="red")
    legend([fz_name(k), "Threshold"], 'Location', 'southeastoutside')
    end
    y_hs_x=[ -350 100;  -1250 200;  -200 2000];
    y_hs_y= [-1700 750 ;  -4000 600; -2500 3000 ];
    y_preppost_x=[-401 250; -150 200; -550 1000];
    y_preppost_y=[-2500 600; -4500 450; -2000 1400];
    figure
    hold on
    for k=1:3
    sgtitle("Post-FES - Forceplate Centre of Preassure x", fontsize=38, FontWeight='bold')
    subplot(3,1,k)    
    plot(tlab, Copx.(k), LineWidth=3)
    grid on
grid minor
    ylabel(["Position of"," the CoP [mm]"], fontsize=18, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
     legend(copx_name(k), 'Location', 'southeastoutside')
        disp(min(Copx.(k))) 
        disp( max(Copx.(k)))
     %xlim([1 6]) %helt
%    xlim([2 13]) %scis
    % ylim([y_hs_x(k,1) y_hs_x(k,2)]) % scis-helt
    ylim([y_preppost_x(k,1) y_preppost_x(k,2)])
     % xlim([2 11]) %prep
    xlim([4 10.5]) %postp
    
    end
    figure
    hold on
    sgtitle("Post-FES - Forceplate Centre of Preassure y", fontsize=38, FontWeight='bold')
    for k=1:3
    subplot(3,1,k)
     plot(tlab, Copy.(k), LineWidth=3)
     grid on
grid minor
     ylabel(["Position of"," the CoP [mm]"], fontsize=18, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
     legend(copy_name(k), 'Location', 'southeastoutside')
    %xlim([1 6]) %helt
%    xlim([2 13]) %Scis
    %ylim([y_hs_y(k,1) y_hs_y(k,2)]) %Scis-helt
    ylim([y_preppost_y(k,1) y_preppost_y(k,2)])
    %xlim([2 11]) %prep
    xlim([4 10.5]) %postp
   disp(min(Copy.(k)) )
       disp(max(Copy.(k)))
    end

    %INITIALIZING THE MASK FOR THE CUT
msk=zeros(size(Fz));
    %CREATING THE MASK
for k=1:3
    th= Fz.(k)<(meanf(k,1)-c(k,1)*stdf(k,1));
    msk(:,k)=th;
end


%CONVERTING THE MASK TO LOGICAL
msk=logical(msk);
%INITIALIZING VECTOR FOR THE LENGTH OF THE MASK
t=zeros(3,1);
%GETTING THE LENGTH OF THE MASK
for k=1:3
tp=length(find(msk(:,k)==1));
t(k,1)=tp;
end

%STANDARDIZATION OF THE MASK
for k=1:3
   b=find(msk(:,k)==1);% find the values where the mask is 1
   beg=b(1); %get the first position where the mask is 1
   egg=b(end)+1; %get the last value and add 1 to get the original length
  if (egg-beg)>min(t) % if the begining - end index is longer than the min length of the shortest mask
      diff=(egg-beg-min(t))/2; %than we calculate the difference and split it half
      if (diff-fix(diff))>0 %the integer part substracted from the whole if we have decimal than this part runs
         diff=fix(diff);
         %msk(1:(beg+diff-1),(egg-diff):length(msk))=0;
         msk(1:beg+diff,k)=0; % distribute the difference between the end and the beginning
         msk(egg-diff:length(msk),k)=0;
      else
        msk(1:(beg+diff),k)=0; %same
        msk((egg-diff+1):length(msk), k)=0;
      end
      
  end
end

% TIMETABLE INITIALIZATION
timetable=zeros(min(t),size(msk,2)); 
%TIMETABLE CREATION
for k=1:3
timetable(:,k)=find(msk(:,k)==1);
end

% DIVISION BY THE SAMPLING FREQUENCY
indextable=timetable; 
timetable=timetable/1200;
% PLOT THE TIMETABLE AND THE MASKED VALUES
tb_EM=zeros(min(t),size(msk,2)*size(EM,2));
tb_fcp=zeros(min(t),size(msk,2)*3);
tb_EM_names=strings(1,size(EM,2)*size(msk,2));

figure
    hold on
    for k=1:3 
        hold on
    subplot(3,1,k)
    hold on
    sgtitle("Post-FES - Cut Force component z", fontsize=38, FontWeight='bold')
    plot( timetable(:,k),Fz.(k)(msk(:,k)),LineWidth=3)
    grid on
grid minor
    ylabel("Froce [N]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    legend(fz_name(k), 'Location', 'southeastoutside')
    xlim([timetable(1,k) timetable(end,k)])
    end

    figure
    hold on
    for k=1:3 
        hold on
    subplot(3,1,k)
    hold on
    sgtitle("Post-FES - Cut Centre of Preassure", fontsize=38, FontWeight='bold')
         plot( timetable(:,k),Copx.(k)(msk(:,k)), LineWidth=3)
         plot( timetable(:,k),Copy.(k)(msk(:,k)), LineWidth=3)
        grid on
grid minor
        ylabel(["Position of"," the CoP [mm]"], fontsize=18, FontWeight='bold')
      
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
    %legend(fz_name(k), 'Location', 'southeast')
     legend([copx_name(k), copy_name(k) ], 'Location', 'southeastoutside')
    xlim([timetable(1,k) timetable(end,k)])
        tb_fcp(:,k)=Fz.(k)(msk(:,k));
        tb_fcp(:,3+k)=Copx.(k)(msk(:,k));
        tb_fcp(:,6+k)=Copy.(k)(msk(:,k));
        for j=1:size(EM,2)
            tb_EM(:,((j-1)*size(EM,2))+k)=EM.(j)(msk(:,k));
            tb_EM_names(:,((j-1)*size(EM,2))+k)=EM.Properties.VariableNames{j}+string(k);
        end
    end
    
    %PLOT THE CUTTED KINEMATICS VALUES

    figure 
    sgtitle("Post-FES - Kinematic tangent cut left",fontsize=38, FontWeight='bold')
    subplot(3,1,1)
    
    plot(timetable(:,1), filtfilt(d1, left_tangent(msk(:,1))), LineWidth=2)
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,1) timetable(end, 1)])
    subplot(3,1,2)
    plot( timetable(:,2), filtfilt(d1, left_tangent(msk(:,2))), LineWidth=2) 
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,2) timetable(end, 2)])
    subplot(3,1,3)
    plot(timetable(:,3), filtfilt(d1, left_tangent(msk(:,3))),LineWidth=2)
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,3) timetable(end, 3)])
    
    avg_cut_left=0;

    figure 
    sgtitle("Post-FES - Kinematic tangent cut right",fontsize=38, FontWeight='bold')
    subplot(3,1,1)
    plot(timetable(:,1), filtfilt(d1, right_tangent(msk(:,1))),LineWidth=2)
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,1) timetable(end, 1)])
    subplot(3,1,2)
    plot(timetable(:,2), filtfilt(d1, right_tangent(msk(:,2))), LineWidth=2) 
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,2) timetable(end, 2)])
    subplot(3,1,3)
    plot(timetable(:,3), filtfilt(d1, right_tangent(msk(:,3))),LineWidth=2)
    grid on
grid minor
    ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,3) timetable(end, 3)])
    avg_cut_right=0;
    % EXCEL GENERATING
    tb_fcp_names=[string(Fz.Properties.VariableNames), string(Copx.Properties.VariableNames), string(Copy.Properties.VariableNames)];
    filename=['cut_data_' f_name];
    dtb_fcp=array2table(tb_fcp);
    dtb_emg=array2table(tb_EM);
    dtb_emg.Properties.VariableNames=tb_EM_names;
    dtb_fcp.Properties.VariableNames=tb_fcp_names;
    for k=0:2

    figure
    if k+1==1
    sgtitle("Post-FES - Cut EMG data of the Vastus Lateralis", fontsize=38, FontWeight='bold')
    elseif k+1==2
     sgtitle("Post-FES - Cut EMG data of the Gastrocnemius", fontsize=38, FontWeight='bold')
    else
        sgtitle("Post-FES - Cut EMG data of the Soleus", fontsize=38, FontWeight='bold')
    end
    for j=1:3
    hold on
    subplot(3,1,j)
    plot(timetable(:,j),dtb_emg.(j+k*3), LineWidth=2)
    grid on
grid minor
    ylabel("sEMG [V]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
    xlim([timetable(1,j) timetable(end, j)])
    ylim([0 1])
    end    
    end
  
    [fcp_tanget, emg_side]=FPA_forceplate(dtb_fcp, EM, timetable);
    [dtb_kin]=kin_cut(upsampled_kin, msk,min(t));
      new_r=[nan(right(1)-1,1); right_tangent];
      new_l= [nan(left(1)-1,1); left_tangent];
    
      if length(new_r)==length(new_l)
      comb= [new_r , new_l];
      else
          nanright=length(upsampled_kin.LHEEx)-length(new_r);
          nanleft=length(upsampled_kin.LHEEx)-length(new_l);
          new_ll=[new_l; nan(nanleft,1)];
          new_rr=[new_r;nan(nanright,1)];
          comb=[new_rr, new_ll];
      end
    dtb_msk=array2table(msk);
    dtb_tan=array2table(comb);
    dtb_time=array2table(timetable);
    dtb_msk.Properties.VariableNames=[ "mask_1", "mask_2","mask_3"];
    dtb_tan.Properties.VariableNames= ["Right_tangent", "Left_tangent"];
    dtb_time.Properties.VariableNames= ["time_1", "time_2", "time_3"];

t1=[timetable(1,1), timetable(end,1), timetable(end,1), timetable(1,1)];
t2=[timetable(1,2), timetable(end,2), timetable(end,2), timetable(1,2)];
t3=[timetable(1,3), timetable(end,3), timetable(end,3), timetable(1,3)];
c1=[0.9725490196078431, 0.8,0.6823529411764706];
c2=[0.9529411764705882,0.6509803921568628, 0.44313725490196076];
c3=[0.9294117647058824, 0.49019607843137253,0.19215686274509805 ];

figure
hold on
plot((1:1:length(new_l))/1200,new_l, LineWidth=2)
xline(timetable(1,1), Color='red', LineWidth=2);
xline(timetable(end,1), Color='red', LineWidth=2);
xline(timetable(1,2), Color='red', LineWidth=2);
xline(timetable(end,2), Color='red', LineWidth=2);
xline(timetable(1,3), Color='red', LineWidth=2);
xline(timetable(end,3), Color='red', LineWidth=2);
lr=ylim;
y1=[lr(1),lr(1), lr(2), lr(2)];
fill(t1,y1,c1, 'FaceAlpha', 0.5);
fill(t2,y1,c2, 'FaceAlpha', 0.5);
fill(t3, y1,c3,'FaceAlpha', 0.5);
grid on
grid minor
title("Post-FES - Full kinematic tangent left", fontsize=38, FontWeight='bold')
ylabel("Degree [°]", fontsize=22, FontWeight='bold')
    xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
   xlim([(left(1,1)-1)/1200 left(end)/1200])
     
figure
hold on
plot((1:1:length(new_r))/1200,new_r, LineWidth=2)
xline(timetable(1,1), Color='red', LineWidth=2);
xline(timetable(end,1), Color='red', LineWidth=2);
xline(timetable(1,2), Color='red', LineWidth=2);
xline(timetable(end,2), Color='red', LineWidth=2);
xline(timetable(1,3), Color='red', LineWidth=2);
xline(timetable(end,3), Color='red', LineWidth=2);
lr2=ylim;
y1=[lr2(1),lr2(1), lr2(2), lr2(2)];
fill(t1,y1,c1, 'FaceAlpha', 0.5);
fill(t2,y1,c2, 'FaceAlpha', 0.5);
fill(t3, y1,c3,'FaceAlpha', 0.5);
grid on
grid minor
title("Post-FES - Full kinematic tangent right", fontsize=38, FontWeight='bold')
ylabel("Degree [°]", fontsize=22, FontWeight='bold')
  xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 22);
   xlim([(right(1,1)-1)/1200 right(end)/1200])

figure
subplot(3,1,1)
plot(dtb_time.time_1, dtb_kin.LHEEz1, LineWidth=2)
title("Post-FES - Left heel z component 1 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
xlabel ("Time [s]", fontsize=22, FontWeight='bold')
set(gca, 'FontSize', 18);
xlim([timetable(1,1) timetable(end, 1)])
grid on
grid minor
subplot(3,1,2)
plot(dtb_time.time_2, dtb_kin.LHEEz2, LineWidth=2)
title("Post-FES - Left heel z component 2 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
xlabel ("Time [s]", fontsize=22, FontWeight='bold')
set(gca, 'FontSize', 18);
xlim([timetable(1,2) timetable(end, 2)])
grid on
grid minor
subplot(3,1,3)
    plot(dtb_time.time_3, dtb_kin.LHEEz3, LineWidth=2)
    grid on
grid minor
title("Post-FES - Left heel z component 3 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
  xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
    xlim([timetable(1,3) timetable(end, 3)])

    figure
subplot(3,1,1)
plot(dtb_time.time_1, dtb_kin.RHEEz1, LineWidth=2)
grid on
grid minor
title("Post-FES - Right heel z component 1 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
  xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
    xlim([timetable(1,1) timetable(end, 1)])
subplot(3,1,2)
    plot(dtb_time.time_2, dtb_kin.RHEEz2, LineWidth=2)
    grid on
grid minor
title("Post-FES - Right heel z component 2 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
  xlabel ("Time [s]", fontsize=18, FontWeight='bold')
    set(gca, 'FontSize', 18);
    xlim([timetable(1,2) timetable(end, 2)])
subplot(3,1,3)
    plot(dtb_time.time_3, dtb_kin.RHEEz3, LineWidth=2)
    grid on
grid minor
title("Post-FES - Right heel z component 3 forceplate", fontsize=38, FontWeight='bold')
ylabel(["Position on"," the z axis [mm]"], fontsize=18, FontWeight='bold')
  xlabel ("Time [s]", fontsize=22, FontWeight='bold')
    set(gca, 'FontSize', 18);
    xlim([timetable(1,3) timetable(end, 3)])

%LINEAR REGRESSION
X=[dtb_emg.r_Vast_Lat_1, dtb_emg.r_Gastrocnem_1, dtb_emg.r_Soleus1];
mdl=fitlm(X, dtb_tan.Right_tangent(msk(:,1)));
    
    writetable(dtb_fcp,filename, 'Sheet','fcp_data')
    writetable(dtb_emg,filename, 'Sheet','emg_data')
    writetable(dtb_kin, filename, 'Sheet', 'kin_data')
    writetable(dtb_tan, filename, 'Sheet', 'tan_data')
    writetable(dtb_time, filename, 'Sheet', 'time_data')
    writetable(dtb_msk, filename, 'Sheet', 'msk_data')

    %% Trash
% CENTRE OF PREASSURE
% dz=82.5;
% dy=400;
% dx=600;
% copx1=(-1*(fcp.fp1_M_y)+fcp.fp1_F_x*dz)./fcp.fp1_F_z;
% copy1=(fcp.fp1_M_x+fcp.fp1_F_y*dz)./fcp.fp1_F_z;
% % copy1=filtfilt(d1, copy1);
% % copx1=filtfilt(d1, copx1);
% Copy.fp1_CoP_y=copy1;
% Copx.fp1_CoP_x=copx1;
% figure
% hold on
% plot(copx1);
% plot(copy1);
