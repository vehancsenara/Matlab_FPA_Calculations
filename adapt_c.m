function [c] = adapt_c(fz, stdfz, meanfz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
len=1:length(fz.fp1_F_z);
c=zeros(3,1);

for k=1:3
    minima=abs(min(fz.(k)));
    locmax_0=islocalmax(fz.(k));
    interaction=fz.(k)(locmax_0)<-400;%only interaction
    logimax=len(locmax_0);
    filtered_idx=logimax(interaction);
    maximas=abs(fz.(k)(filtered_idx));
    diff=minima-maximas;
    max_diff=max(diff);
    norm_diff=abs((diff/max_diff)-0.5);
    min_diff=min(norm_diff);
    x=(minima-diff(norm_diff==min_diff))*0.85;
    c(k,1)=((x)-meanfz(k,1))/stdfz(k,1);
end

end