function [ordered_fz, ordered_copy,ordered_copx] = order_fz(fz, copy, copx, freq)
ordered_fz=array2table(zeros(size(fz)));
ordered_copy=array2table(zeros(size(copy)));
ordered_copx=array2table(zeros(size(copx)));
ordered_copx.Properties.VariableNames=copx.Properties.VariableNames;
ordered_copy.Properties.VariableNames=copy.Properties.VariableNames;
ordered_fz.Properties.VariableNames=fz.Properties.VariableNames;
fmins=zeros(width(fz),2);
for k=1:width(fz)
    [~, minidx]=min(fz.(k));
    fmins(k, 1)=minidx;
    fmins(k,2)=k;
end
sorted_fmins=sortrows(fmins,1);
for k=1:width(fz)
    ordered_fz.(k)=fz.(sorted_fmins(k,2));
    ordered_copx.(k)=copx.(sorted_fmins(k,2));
    ordered_copy.(k)=copy.(sorted_fmins(k,2));
end
if freq ~=length(ordered_fz.(1))
    upsampled_fz=zeros(freq, width(ordered_fz));
    upsampled_copx=zeros(freq, width(copx));
    upsampled_copy=zeros(freq, width(copy));
for k=1:width(ordered_fz)
    len_fz=1:1:length(ordered_fz.(k));
    orig_fz=ordered_fz.(k);
    upsampletime_fz=linspace(len_fz(1), len_fz(end), length(len_fz)*(freq/length(len_fz)));
    upsampled_fz(:,k)=interp1(len_fz, orig_fz, upsampletime_fz, 'linear');

    len_copx=1:1:length(ordered_copx.(k));
    orig_copx=ordered_copx.(k);
    upsampletime_copx=linspace(len_copx(1), len_copx(end), length(len_copx)*(freq/length(len_copx)));
    upsampled_copx(:,k)=interp1(len_copx, orig_copx, upsampletime_copx, 'linear');

    len_copy=1:1:length(ordered_copy.(k));
    orig_copy=ordered_copy.(k);
    upsampletime_copy=linspace(len_copy(1), len_copy(end), length(len_copy)*(freq/length(len_copy)));
    upsampled_copy(:,k)=interp1(len_copy, orig_copy, upsampletime_copy, 'linear');
end
ordered_fz=array2table(upsampled_fz);
ordered_fz.Properties.VariableNames=fz.Properties.VariableNames;

ordered_copx=array2table(upsampled_copx);
ordered_copx.Properties.VariableNames=copx.Properties.VariableNames;

ordered_copy=array2table(upsampled_copy);
ordered_copy.Properties.VariableNames=copy.Properties.VariableNames;
end

end