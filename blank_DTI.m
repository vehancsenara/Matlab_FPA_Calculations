clear;
%beolvasás

data = niftiread("4Ddwi_b1000.nii"); %gondolom innet jon az ezer
mask = niftiread("brain_mask.nii");
grads = readmatrix("grad_dirs.txt"); %dlmread-re cserélhető

%pre-alloc
num_scan = length(grads);
scan_0 = zeros(size(data(:,:,:,1)));
num_0 = 0;
scan = zeros([size(scan_0),num_scan]);
gr = zeros(size(grads(1,:)));
i=0;

%adatok kettéválasztása
for k=1:num_scan
    if(isequal(grads(k,:),[0,0,0])) %ha a gradiens 0
        num_0 = num_0 + 1;
        scan_0 = scan_0 + data(:,:,:,k);
    else
        i=i+1;
        scan(:,:,:,i)=data(:,:,:,k); %4D adat tárolása
        gr(i,:)=grads(k,:); %hozzá tartozó gradiens
    end
end
scan_0 = scan_0/num_0;
clear data grads
%pre-alloc
scan_s = size(scan_0);
b = zeros([3,3,size(gr,1)]);
scan_ab = zeros([scan_s,2]);
DifT = zeros([scan_s,6]);
Y = zeros([scan_s,3]);
FA = zeros(scan_s);
ADC = zeros(scan_s); 
V_dir = zeros([scan_s,3]);

for k = 1:size(gr,1)
    b(:,:,k)=gr(k,:)'*gr(k,:)*1000;
    scan_ab(:,:,:,k)=log((scan(:,:,:,k)./scan_0)+eps);
end

b_vec = squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

th_wm = 0.1; %fehérállomány threshold

for x = 1:scan_s(1)
    for y = 1:scan_s(2)
        for z = 1:scan_s(3)
            if(mask(x,y,z))

                Z=squeeze(scan_ab(x,y,z,:));
                M=b_vec\Z;

                DifTensor=[M(1) M(2) M(3);M(2) M(4) M(5); M(3) M(5) M(6)];

                [SajV,D]=eig(DifTensor);
                SajE=diag(D);

                [~,indx]=sort(SajE);
                SajE=SajE(indx);
                SajV=SajV(:,indx);
                Saje_0=SajE;

                if((SajE(1)<0)&&(SajE(2)<0)&&(SajE(3)<0))
                    SajE=abs(SajE);
                end

                if (SajE(1)<=0)
                    SajE(1)=eps;
                end

                if (SajE(2)<=0)
                    SajE(2)=eps;
                end

                if (SajE(3)<=0)
                    SajE(3)=eps;
                end

                ADC_V=(sum(SajE, 'all'))/3;
                
                szaml=sqrt((SajE(1)-SajE(2)).^2 + (SajE(1)-SajE(3)).^2 +(SajE(2)-SajE(3)).^2);
                nevz= sqrt(SajE(1).^2+SajE(2).^2+SajE(3).^2);
                FA_V=(1/sqrt(2))*(szaml./nevz);

                ADC(x,y,z)=ADC_V;
                Y(x,y,z,:)=SajE;
                DifT(x,y,z,:)=[DifTensor(1:3), DifTensor(5:6), DifTensor(9)];
                if (FA_V>th_wm)
                    FA(x,y,z)=FA_V;
                    V_dir(x,y,z,:)=SajV(:,end)*Saje_0(end);
                end

            end

        end
    end
end