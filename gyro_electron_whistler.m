
%%
path='E:\projects\Elec_chorus\Elec_chorus\FSD.txt';
[t,vx,vy,vz,enem,Bx,By,z]=textread(path,'%f %f %f %f %f %f %f %f','delimiter',' ');

E1=1/2*9.1*1e-31*(vx.^2+vy.^2+vz.^2)/1.6*1e19*1e10;
vxpla=0;
vypla=0;
vzpla=0;
vz=vz+vzpla/1e5;
vx=vx+vxpla/1e5;
vy=vy+vypla/1e5;
E=1/2*9.1*1e-31*(vx.^2+vy.^2+vz.^2)/1.6*1e19*1e10;

%% load data part 2
maxindextotal=[];
maxindextotal(1)=1;
for i=2:length(t)
    if abs(t(i)-t(i-1))>0.05
        maxindextotal=[maxindextotal,i];
        else if (t(i)-t(i-1)) == 0
                maxindextotal=[maxindextotal,i];
            else if t(i)>t(i-1)
                maxindextotal=[maxindextotal,i];
                end
            end
    end
end
[Vphi,Vperp]=cart2pol(vx,vy);
[Bphi,Bperp]=cart2pol(Bx,By);
Vphi=Vphi*180/pi;
p=find(Vphi<0);
Vphi(p)=Vphi(p)+360;
Bphi=Bphi*180/pi;
p=find(Bphi<0);
Bphi(p)=Bphi(p)+360;
theta=Vphi-Bphi;
poshigh=find(theta>=180);
theta(poshigh)=theta(poshigh)-360;
neglow=find(theta<=-180);
theta(neglow)=theta(neglow)+360;
%%
for ener=1
% timepoint=36;
timepoint=300;
gyro=36;
for pt=1
xxx=1*timepoint*gyro*(ener-1);
maxindex=maxindextotal(gyro*timepoint*(pt-1)+1+xxx:gyro*timepoint*pt+xxx); % the end channel no +1
datasave50eV=cell(gyro,timepoint);
for i=1:timepoint*gyro 
%     if i==timepoint*gyro && pt==1 && ener==32
%     if i==timepoint*gyro && pt==1 && ener==12
    if i==timepoint*gyro && pt==1 && ener==1
        index=maxindex(i):length(t);
    else 
        index=maxindex(i):maxindex(i+1)-1;
    end
    y=ceil(i/gyro);
    x=mod(i,gyro);
    if x==0
        x=gyro;
    end
%       datasave50eV{x,y}={-t(index(end)),E(index(end)),theta(index(end)),Vperp(index(end))*1e5,vx(index(end))*1e5,vy(index(end))*1e5,vz(index(end))*1e5,Bx(index(end)),By(index(end)),E1(index(end))};
    datasave50eV{x,y}={-t(index),E(index),theta(index),Vperp(index)*1e5,vx(index)*1e5,vy(index)*1e5,vz(index)*1e5,Bx(index),By(index),E1(index),z(index)};
end
% c_eval('save(''datasave30eV_PA130_gyro30_360_amp=1.5nT_B0_2nT_f0.45_5phase_vphae400_real.mat'',''datasave50eV'')',2.5*pt)
% c_eval('save(''datasave17eV_PA40_gyro30_360_amp=1.5nT_B0_2nT_f0.45_5phase_vphae400_tp41.mat'',''datasave50eV'')',2.5*pt)
c_eval('save(''datasave39eV_PA140_gyro36_amp=1.5nT_B0_2nT_fvary_vphae400_2_bochang_increase_tp300_modify_gaussian_phasechanged.mat'',''datasave50eV'')',10*pt+10)
end
end
%%  calculate 3-D distribution
psd_=[];
psd=[];
dataplot=datasave50eV;
vxpla=0e3;
vypla=0e3;
vzpla=450e3;
mp=9.1*1e-31;
qi=1.6*1e-19;
Tperp=10;
Tpara=20;
const=1.8*1e6*(mp/2/pi/Tperp/qi)^(2/2)*(mp/2/pi/Tpara/qi)^(1/2)*1e13;
f=@(vx,vy,vz)exp(-mp/2*(vx-vxpla).^2/Tperp/qi-mp/2*(vy-vypla).^2/Tperp/qi-mp/2*(vz-vzpla).^2/Tpara/qi);
for i=1:gyro
    for j=1:timepoint
        vx_=cell2mat(dataplot{i,j}(5)); 
        vx_=vx_(end);
        vy_=cell2mat(dataplot{i,j}(6));
        vy_=vy_(end);
        vz_=cell2mat(dataplot{i,j}(7));
        vz_=vz_(end);
       psd(i,j)=f(vx_,vy_,vz_);
    end
end
psd_=[psd;zeros(1,timepoint)];
psd_=[psd_,zeros(gyro+1,1)];
h=pcolor(log10(psd_*const));
set(h,'LineStyle','none');

epoch = EpochUnix(iso2epoch('2002-03-04T09:30:00Z')+(0:0.06:299*0.06));
load('B_modify_gaussian.mat');
Bz=zeros(length(Bx_),1);
Bgse=TSeries(epoch,[Bx_,By_,Bz])
Bphi=transphase2(Bgse,1);
hold on
yyaxis right
plot(Bphi)
ylim([0,360])
colormap('jet')
