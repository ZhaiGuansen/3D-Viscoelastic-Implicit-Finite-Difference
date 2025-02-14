% Cylindrical beamformer  transform method
% 
function [Spec,f,vel]=BeamformerFV(seis)
Ns=1;
[nt,nx]=size(seis.data);
dt=seis.dt;
t = 0:dt:dt*(nt-1);
h=seis.h;

flow=seis.flow;
fhigh=seis.fhigh;
vlow=seis.vlow;
vhigh=seis.vhigh;
dv=seis.dv;
vel=vlow:dv:vhigh;
vlen=length(vel);

data=seis.data;

%[nt,nq] = size(m);
nh = length(h);

nfft = 4*(2^nextpow2(nt));

M = fft(data,nfft,1);


D = zeros(nfft,nh);

i = sqrt(-1);

ilow  = floor(flow*dt*nfft)+1;
if ilow < 1; ilow=1; end;
ihigh = floor(fhigh*dt*nfft)+1;
if ihigh > floor(nfft/2)+1; ihigh = floor(nfft/2)+1; end

D = zeros(nh,nh);
for ifreq=ilow:1:ihigh
    x = M(ifreq,:)';
    for j=1:nh
       for jj=1:nh
          R(j,jj)=x(j)*conj(x(jj));
          R(jj,j)=x(jj)*conj(x(j));
       end
    end
   
      [V1,S]=eig(R);
      Rinv=zeros(nh,nh);
      for jm=1:nh-0
        V2=V1(:,jm);
        Rinv=Rinv+V2*V2';%.*S(jm,jm);
      end
        
     % if(ifreq==30) figure; imagesc(abs(R)); figure;imagesc(abs(Rinv));disp( S);end

      for iv=1:vlen
       v= vel(iv);
       f0 = 2.*pi*(ifreq-1)/nfft/dt;
       for ix=1:nh
          L(ix) = exp(-i*f0*(h(ix)/v));  %plane wave
       y0=besselh(0,f0*h(ix)/v+0.00000001);   %  Cylindrical beamformer kr Hankel fuction
       L(ix) = exp(-i*(phase(y0)+pi./4));  %
       end
             
       PP=(L)*(R)*((L)') ;      %  beamformer  transform
   %   PP=1./((L)*(Rinv)*((L)')) ; %FV-MUSIC
       P(ifreq-ilow+1,iv)=abs(PP);
    end

    
   
  f1(ifreq-ilow+1)=f0/2/pi;  
 
 
 end
 
% figure
%  imagesc(vel,f1,abs(P));
% 
  f=f1;
  Spec=P;

end

