     path=pwd;
     name1='all3';
     name2='\output';
     num='5';
     d='x';
     name3='.dat';

     file=strcat(path,'\',name1,name2,num,'_',d,name3);

     data=readdat(file,800);

     dt = 0.1;%单位ms
     dx = 1;
     
     data=data./max(max(data));
     data=data(:,70:130);
     [nt,nx]=size(data);
     x=1*dx:dx:nx*dx;
     t=1*dt:dt:nt*dt;
     wigb(data,1,x,t);
     
     
     minOff=dx*40;%炮点
     h = minOff:dx:(nx-1)*dx+minOff;%偏移距
     h=abs(h);
     flow =10;fhigh =350;
     vlow=500;vhigh=4000;dv=10;
     vel = vlow:dv:vhigh;
     
     
     seis.data=data;
     seis.dt=dt*1e-3;
     seis.nt=nt;
     seis.dx=dx;
     seis.nx=nx;
     seis.minOff=minOff;
     seis.flow=flow;
     seis.fhigh=fhigh;
     seis.vlow=vlow;
     seis.vhigh=vhigh;
     seis.dv=dv;
     seis.h=h;
     seis.vel=vel;
     
 % beamformer
     [spec,f,k]=BeamformerFV(seis);
  figure;
  pcolor(f,vel,abs(spec)');shading interp;
  colorbar;
  %clim([0,4000])
  xlabel('Frequency/Hz'); ylabel('Velocity/m/s')



% %谱增强

  D=Spec_Enhence(abs(spec),10,20,150,1);
  figure;
  pcolor(f,vel,abs(D)');shading interp;
  colorbar;
  xlabel('Frequency/Hz'); ylabel('Velocity/m/s')
%    
% figure
% imagesc(k,f,D);
% %set(gca,'YDir','normal');
% set(gca,'XAxisLocation','top'); 


