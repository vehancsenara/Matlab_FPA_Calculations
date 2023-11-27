function [left_FPA,right_FPA, right_tangent, left_tangent, left, right, upsampled_kin] = FPA_kin(kin, freq)
%azért mert különböző helyeken van nan érték ezért megnézzük azokat az
%indexeket, ahol nem nan az érték, majd veszsük a metszetét és abból
%számoljuk a többit a kivonom vektor tangens hármas
upsampled_kin=zeros(freq, width(kin));
for k=1:width(kin)
    nan_idx=isnan(kin.(k));
    x=1:1:length(kin.(k));
    x_query=x(nan_idx);
    kin.(k)(nan_idx) = interp1(x(~nan_idx), kin.(k)(~nan_idx), x_query, 'linear');
    kinlen=1:1:length(kin.(k));
    orig_kin=kin.(k);
    upsampletime=linspace(kinlen(1), kinlen(end), length(kinlen)*(freq/length(kinlen)));
    upsampled_kin(:,k)=interp1(kinlen, orig_kin, upsampletime, 'linear');
end
upsampled_kin=array2table(upsampled_kin);
upsampled_kin.Properties.VariableNames=kin.Properties.VariableNames;

lefttoe=find(~isnan(upsampled_kin.LTOEx));
leftheel=find(~isnan(upsampled_kin.LHEEx));
left=intersect(lefttoe, leftheel);



leftnx= upsampled_kin.LTOEx(left)-upsampled_kin.LHEEx(left);
leftny= upsampled_kin.LTOEy(left)-upsampled_kin.LHEEy(left);
left_tangent=atand(leftnx./leftny);
left_FPA=mean(left_tangent);

righttoe=find(~isnan(upsampled_kin.RTOEx));
rightheel=find(~isnan(upsampled_kin.RHEEx));
right=intersect(righttoe, rightheel);



rightnx= upsampled_kin.RTOEx(right)-upsampled_kin.RHEEx(right);
rightny= upsampled_kin.RTOEy(right)-upsampled_kin.RHEEy(right);
right_tangent=atand(rightnx./rightny);
right_FPA=mean(right_tangent);


end