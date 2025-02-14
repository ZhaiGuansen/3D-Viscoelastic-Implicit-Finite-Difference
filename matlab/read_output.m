function [size,value_x,value_y,value_z]=read_output(T)
%num=input('选择非固定参数：(x,y,z,t)');
%n=input('输入剩余三个固定参数：');

path=pwd;
name1='\OUTPUT\OUTPUT';
name2='.DAT';
for t=1:T
    fname=strcat(path,name1,num2str(t),name2);
    fid=fopen(fname,"r");
    %line=fgetl(fid);
    size=fscanf(fid,'%d',3);
    if t==1
        value_x=zeros(size(1),size(2),size(3),T);
        value_y=zeros(size(1),size(2),size(3),T);
        value_z=zeros(size(1),size(2),size(3),T);
    end
    % data=fscanf(fid,'%lf',size(1)*size(2)*size(3));
    % value_x(:,:,:,t)=reshape(data,size(1),size(2),size(3));
    % 
    % data=fscanf(fid,'%lf',size(1)*size(2)*size(3));
    % value_y(:,:,:,t)=reshape(data,size(1),size(2),size(3));
    % 
    % data=fscanf(fid,'%lf',size(1)*size(2)*size(3));
    % value_z(:,:,:,t)=reshape(data,size(1),size(2),size(3));


    for z=1:size(3)
        for y=1:size(2)
            value_x(:,y,z,t)=fscanf(fid,'%lf',size(1));
            % for x=1:size(1)
            %     value_x(x,y,z,t)=fscanf(fid,'%lf',1);
            % end
        end
    end
    for z=1:size(3)
        for y=1:size(2)
            value_y(:,y,z,t)=fscanf(fid,'%lf',size(1));
            % for x=1:size(1)
            %     value_y(x,y,z,t)=fscanf(fid,'%lf',1);
            % end
        end
    end
    for z=1:size(3)
        for y=1:size(2)
            value_z(:,y,z,t)=fscanf(fid,'%lf',size(1));
            % for x=1:size(1)
            %     value_z(x,y,z,t)=fscanf(fid,'%lf',1);
            % end
        end
    end
    if mod(t,20)==0
    fprintf('%d/%d\n',t,T);
    end
end
fprintf('读取完成')
end
