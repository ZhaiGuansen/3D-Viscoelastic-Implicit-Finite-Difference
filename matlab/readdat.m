function D = readdat(file,t1)
f=fopen(file,'r');
C=textscan(f, '%f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
data=C{1};
nx=data(1);
nt=data(2);
%t1=1200;
%t1=nt;
i=3;
D=zeros(t1,nx);
for x=1:nx
    for t=1:nt
        if t<=t1
            D(t,x)=data(i);
        end
        i=i+1;
    end
end


end

