clear;
close all;

fcp=readtable('h06_r_1stp1fcp_slow_1.xlsx','Sheet','fcp');
hasMatchfz = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'F_z', 'once'));
hasMatchcx = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'CoP_x', 'once'));
hasMatchcy = ~cellfun('isempty', regexp(fcp.Properties.VariableNames, 'CoP_y', 'once'));
Fz=fcp(:, fcp.Properties.VariableNames(hasMatchfz));
Copx=fcp(:, fcp.Properties.VariableNames(hasMatchcx));
Copy=fcp(:, fcp.Properties.VariableNames(hasMatchcy));

d1=designfilt("lowpassiir",FilterOrder=6, ...
    HalfPowerFrequency=0.1,DesignMethod="butter");

for k=1:3
fz=filtfilt(d1, Fz.(k));
cx=filtfilt(d1, Copx.(k));
cy=filtfilt(d1, Copy.(k));
Fz.(k)=fz;
Copx.(k)=cx;
Copy.(k)=cy;
end


v1=fcp.fp1_F_z;
v2= fcp.fp2_F_z;
v3=fcp.fp3_F_z;
cx1=fcp.fp1_CoP_x;
cy1=fcp.fp1_CoP_y;
cx2=fcp.fp2_CoP_x;
cy2=fcp.fp2_CoP_y;
cx3=fcp.fp3_CoP_x;
cy3=fcp.fp3_CoP_y;


fv1=filtfilt(d1, v1);
fv2=filtfilt(d1, v2);
fv3=filtfilt(d1, v3);
figure
subplot(3,1,1)
hold on
plot(v1, 'LineWidth',2);
plot(fv1, LineWidth=2, LineStyle="--");
subplot(3,1,2)
hold on
plot(v2, 'LineWidth',2);
plot(fv2, LineWidth=2, LineStyle="--");
subplot(3,1,3)
hold on
plot(v3, 'LineWidth',2);
plot(fv3, LineWidth=2, LineStyle="--");
%mean-std 
fvm1=mean(fv1(1:500));
fvs1=std(fv1(1:500));
fvm2=mean(fv2(1:500));
fvs2=std(fv2(1:500));
fvm3=mean(fv3(1:500));
fvs3=std(fv3(1:500));
%histograms
his1=fv1(1:500);
his2=fv2(1:500);
his3=fv3(1:500);
%figures HISTOGRAM
figure
sgtitle("Distribution of the Noise of the Force plate",fontsize=34, FontWeight='bold')
subplot(3,1,1)
histfit(his1);
ylabel("Frequency", fontsize=18, FontWeight='bold')
    xlabel ("Noise sample [N]", fontsize=18, FontWeight='bold')
    set(gca, 'FontSize', 14);
subplot(3,1,2)
histfit(his2);
ylabel("Frequency", fontsize=18, FontWeight='bold')
    xlabel ("Noise sample [N]", fontsize=18, FontWeight='bold')
    set(gca, 'FontSize', 14);
subplot(3,1,3)
histfit(his3)
ylabel("Frequency", fontsize=18, FontWeight='bold')
    xlabel ("Noise sample [N]", fontsize=18, FontWeight='bold')
    set(gca, 'FontSize', 14);
pd=fitdist(his2, 'Normal');

trashhold1=fvm1-250*fvs1;
trashhold2=fvm2-250*fvs2;
trashhold3=fvm3-250*fvs3;
figure
subplot(3,1,1)
hold on
yline(trashhold1, LineWidth=1, Color="red")
plot(fv1)
subplot(3,1,2)
hold on
yline(trashhold2,LineWidth=1, Color="red")
plot(fv2)
subplot(3,1,3)
hold on
yline(trashhold3, LineWidth=1, Color="red")
plot(fv3)

cutcx1=zeros(size(fv1));
cutcy1=zeros(size(fv1));
cutcx2=zeros(size(fv2));
cutcy2=zeros(size(fv2));
cutcx3=zeros(size(fv3));
cutcy3=zeros(size(fv3));
msk1=fv1<trashhold1;
msk2=fv2<trashhold2;
msk3=fv3<trashhold3;
cutcx1(msk1)=filtfilt(d1,cx1(msk1));
cutcy1(msk1)=filtfilt(d1,cy1(msk1));
cutcx2(msk2)=filtfilt(d1,cx2(msk2));
cutcy2(msk2)=filtfilt(d1,cy2(msk2));
cutcx3(msk3)=filtfilt(d1,cx3(msk3));
cutcy3(msk3)=filtfilt(d1,cy3(msk3));
cutter1=zeros(size(fv1));
cutter1(fv1<trashhold1)=fv1(fv1<trashhold1);
cutter2=zeros(size(fv2));
cutter2(fv2<trashhold2)=fv2(fv2<trashhold2);
cutter3=zeros(size(fv3));
cutter3(fv3<trashhold3)=fv3(fv3<trashhold3);


figure
hold on
subplot(3,1,1)
hold on
plot(cutter1)
plot(cutcx1)
plot(cutcy1)
subplot(3,1,2)
hold on
plot(cutter2)
plot(cutcx2)
plot(cutcy2)
subplot(3,1,3)
hold on
plot(cutter3)
plot(cutcx3)
plot(cutcy3)
