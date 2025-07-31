function [size,value_x,value_y,value_z]=read_output(T)

folder = uigetdir('','Please select the save location for the OUTPUT file.');

if folder == 0
    error('The user cancelled the selection.');
end
name1='/OUTPUT';
name2='.DAT';

for t=1:T
    fname=strcat(folder,name1,num2str(t),name2);
    fid=fopen(fname,"r");
    %line=fgetl(fid);
    size=fscanf(fid,'%d',3);
    if t==1
        value_x=zeros(size(1),size(2),size(3),T);
        value_y=zeros(size(1),size(2),size(3),T);
        value_z=zeros(size(1),size(2),size(3),T);
    end


    for z=1:size(3)
        for y=1:size(2)
            value_x(:,y,z,t)=fscanf(fid,'%lf',size(1));
        end
    end
    for z=1:size(3)
        for y=1:size(2)
            value_y(:,y,z,t)=fscanf(fid,'%lf',size(1));
        end
    end
    for z=1:size(3)
        for y=1:size(2)
            value_z(:,y,z,t)=fscanf(fid,'%lf',size(1));
        end
    end
    if mod(t,20)==0
    fprintf('%d/%d\n',t,T);
    end
end
fprintf('finish')
end
