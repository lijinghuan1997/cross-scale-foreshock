path='C:\Users\dell\Desktop\shocklet_code\Console1\Console1\FSD.txt';
[t,vx,vy,vz,enem,Bx,By]=textread(path,'%f %f %f %f %f %f %f');
vxpla=0; %-22e3;
vypla=0;  %11e3;
vzpla=0;  %35e3;
vz=vz+vzpla/1e5;
vx=vx+vxpla/1e5;
vy=vy+vypla/1e5;
E1=1/2*1.67*1e-27*(vx.^2+vy.^2+vz.^2)/1.6*1e19*1e10;
maxindextotal=[];
maxindextotal(1)=1;
for i=2:length(t)
    if t(i)>t(i-1)
        maxindextotal=[maxindextotal,i];
    else if t(i)==max(t)
            maxindextotal=[maxindextotal,i];
        else if (t(i)-t(i-1))<-0.1
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


%% max  2*pi
poshigh=find(theta>=180);
theta(poshigh)=theta(poshigh)-360;
neglow=find(theta<=-180);
theta(neglow)=theta(neglow)+360;

%% 
for ener=1
for pt=1
timepoint=151;
gyro=12;
xxx=timepoint*gyro*(ener-1);
maxindex=maxindextotal(gyro*timepoint*(pt-1)+1+xxx:gyro*timepoint*pt+xxx); % the end channel +1

datasave50eV=cell(gyro,timepoint);
for i=1:timepoint*gyro 
    if i==timepoint*gyro && ener==1
        index=maxindex(i):length(t);
    else 
        index=maxindex(i):maxindex(i+1)-1;
    end
    y=ceil(i/gyro);
    x=mod(i,gyro);
    if x==0
        x=gyro;
    end
   datasave50eV{x,y}={-t(index),E1(index),theta(index),Vperp(index)*1e5,vx(index)*1e5,vy(index)*1e5,vz(index)*1e5,Bx(index),By(index),E1(index)};
end

c_eval('save(''datasave3keV_PA_80_Bw_2nT_B0_3nT_tp_151_WNA_0_GYRO_12_2lambda_vphase450.mat'',''datasave50eV'');',ener*20+20)
end
end
%% gyroplot
Tint=irf.tint('2021-01-07T08:16:00.00Z/2021-01-07T08:17:00.00Z');
c_eval('ePDist_i=mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',1);
index=ePDist_i.data==0;
ePDist_i.data(index)=NaN;
data=squeeze(irf.nanmean(ePDist_i.data,1));
psd=zeros(gyro,timepoint);
mp=1.67*1e-27;
qi=1.6*1e-19;
for ploti=1
        dataplot=datasave50eV;
for i=1:gyro
    for j=1:timepoint
        vx_=cell2mat(dataplot{i,j}(5)); 
        vx_=vx_(end);
        vy_=cell2mat(dataplot{i,j}(6));
        vy_=vy_(end);
        vz_=cell2mat(dataplot{i,j}(7));
        vz_=vz_(end);
        
        v=-[vx_,vy_,vz_]*[0 0 -1; 0 -1 0; -1 0 0];
        [vphi,vperp]=cart2pol(v(1),v(2));
        vphi=vphi*180/pi;
        index=vphi<0;
        vphi(index)=vphi(index)+360;
        pitch=atan2d(vperp,vz_);
        
        phiID=floor(vphi/11.25)+1;
        pitchID=floor(pitch/11.25)+1;
        
        energy=1/2*mp*(v(1)^2+v(2)^2+v(3)^2)/qi;
        energyID=floor(log10(energy/210)/log10(1.1276))+1;
        if energyID<1 | energyID>32
            psd(i,j)=NaN;
        else
            psd(i,j)=data(energyID,phiID,pitchID);
        end
    end
end
psd_=[psd;zeros(1,timepoint)];
psd_=[psd_,zeros(gyro+1,1)];

h=pcolor(log10(psd_))
set(h,'LineStyle','none');
colormap('jet')
% 
time=0:0.01:450;
Bphi_=360-2*pi/45*time/pi*180;
p=find(Bphi_<0);
while ~isempty(p)
    p=find(Bphi_<0);
    Bphi_(p)=Bphi_(p)+360;
    p=find(Bphi_<0);
end
q=find(Bphi_>357);
Bphi_(q)=NaN;
n=find(Bphi_<3);
Bphi_(n)=NaN;


yyaxis right
% % % plot(1.5:0.01*45/13.5:46.5,Bphi_,'linewidth',2)
plot(1.5:0.01*150/450:151.5,Bphi_,'w','linewidth',2)
ylim([0,360])
end



