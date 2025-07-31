function  output_plot3(T,num1,n1,num2,n2,size,value_x,value_y,value_z,cp0,namenum)
%Fix two parameters among x, y, z, and t to obtain a two-dimensional seismic record.
name=['x','y','z','t'];
if num1>num2
    t=num1;
    num1=num2;
    num2=t;
    t=n1;
    n1=n2;
    n2=t;
end
name1=name(num1);
name2=name(num2);
if num1==1
    if num2==2
        namex=name(3);
        namey=name(4);
        l1=size(3);
        l2=T;
        x=reshape(value_x(n1,n2,:,:),l1,l2);
        y=reshape(value_y(n1,n2,:,:),l1,l2);
        z=reshape(value_z(n1,n2,:,:),l1,l2);

    elseif num2==3
        namex=name(2);
        namey=name(4);
        l1=size(2);
        l2=T;
        x=reshape(value_x(n1,:,n2,:),l1,l2);
        y=reshape(value_y(n1,:,n2,:),l1,l2);
        z=reshape(value_z(n1,:,n2,:),l1,l2);

    else
        namex=name(2);
        namey=name(3);
        l1=size(2);
        l2=size(3);
        x=reshape(value_x(n1,:,:,n2),l1,l2);
        y=reshape(value_y(n1,:,:,n2),l1,l2);
        z=reshape(value_z(n1,:,:,n2),l1,l2);

    end


elseif num1==2
    namex=name(1);
    if num2==3
        namey=name(4);
        l1=size(1);
        l2=T;
        x=reshape(value_x(:,n1,n2,:),l1,l2);
        y=reshape(value_y(:,n1,n2,:),l1,l2);
        z=reshape(value_z(:,n1,n2,:),l1,l2);
    else
        namey=name(3);
        l1=size(1);
        l2=size(3);
        x=reshape(value_x(:,n1,:,n2),l1,l2);
        y=reshape(value_y(:,n1,:,n2),l1,l2);
        z=reshape(value_z(:,n1,:,n2),l1,l2);
    end

elseif num1==3
    namex=name(1);
    namey=name(2);
    l1=size(1);
    l2=size(2);
    x=reshape(value_x(:,:,n1,n2),l1,l2);
    y=reshape(value_y(:,:,n1,n2),l1,l2);
    z=reshape(value_z(:,:,n1,n2),l1,l2);


end
a=1:l1;
b=1:l2;
[B,A]=meshgrid(b,a);
if num2~=4
    figure;
    name0=append(name1,'=',num2str(n1),',',name2,'=',num2str(n2));
    set(gcf,'Name',name0);
end
    

    cp1="/output"+namenum;


if num2~=4
    subplot(1,3,1);
    plot3(A,B,x);
    title('x');xlabel(namex);ylabel(namey);
end
    cp=cp0+cp1+"_x.dat";
    fid=fopen(cp,"w");
    fprintf(fid,"%d %d\n",l1,l2);
    for i=1:l1
        for j=1:l2
            fprintf(fid,"%g ",x(i,j));
        end
        fprintf(fid,"\n");
    end
    fclose(fid);
if num2~=4
    subplot(1,3,2);
    plot3(A,B,y);
    title('y');xlabel(namex);ylabel(namey);
end
    cp=cp0+cp1+"_y.dat";
    fid=fopen(cp,"w");
    fprintf(fid,"%d %d\n",l1,l2);
    for i=1:l1
        for j=1:l2
            fprintf(fid,"%g ",y(i,j));
        end
        fprintf(fid,"\n");
    end
    fclose(fid);
if num2~=4
    subplot(1,3,3);
    plot3(A,B,z);
    title('z');xlabel(namex);ylabel(namey);
    set(gcf,'Position',[50 480 1500 500])
end
    cp=cp0+cp1+"_z.dat";
    fid=fopen(cp,"w");
    fprintf(fid,"%d %d\n",l1,l2);
    for i=1:l1
        for j=1:l2
            fprintf(fid,"%g ",z(i,j));
        end
        fprintf(fid,"\n");
    end
    fclose(fid);

end


