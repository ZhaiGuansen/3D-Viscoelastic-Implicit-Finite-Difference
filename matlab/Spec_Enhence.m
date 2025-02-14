%Dispersion Spectrum enhencement processing 
% data, f-k or f-v amplitude spectrum  data
% D-- enhanceing Spectrm f-v spectrum 
% m(mf)- Frequncy numbers
% n(nv) - velocity or wavenumbers;
% a=2,b=20   variable alpha for each frequecy
function D=Spec_Enhence(data,fmin,a,b,span)
[m,n]=size(data);
for i=1:m
    alpha(i)=a+b/(i+fmin);
end

for i=1:m
 for j=1:n
     kdata(j)=data(i,j);
 end
 sdata=smooth(kdata,span);
 cmax(i)=max(sdata);
  for j=1:n
      C(i,j)=data(i,j)./(cmax(i)+0.000001);
      data_en(i,j)=C(i,j).^alpha(i);
   %   data_en(i,j)=C(i,j).*data(i,j);
  end
  
end

D=data_en;

end
