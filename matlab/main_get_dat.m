
nt=200;% number of time steps
%Survey line position parameters
x0=22;
y0=22;
z0=22;
 [size,value_x,value_y,value_z]=read_output(nt);
cp0=pwd;
name="matlab_output";
cp0=fullfile(cp0,name);
if ~exist(cp0,'dir') 
    mkdir(cp0)
end

  output_plot3(nt,1,x0,2,y0,size,value_x,value_y,value_z,cp0,"1");
  output_plot3(nt,1,x0,3,z0,size,value_x,value_y,value_z,cp0,"2");
  output_plot3(nt,2,y0,3,z0,size,value_x,value_y,value_z,cp0,"3");
  %output_plot3(nt,1,90,3,z0,size,value_x,value_y,value_z,cp0,"4");
  %output_plot3(nt,2,y0,3,z0+6,size,value_x,value_y,value_z,cp0,"5");
x1=2;
x2=y0;
for t=40:20:200
    cp=cp0+"/timeshot/"+num2str(x1)+","+num2str(x2)+"/"+num2str(t);
    if ~exist(cp,'dir') 
        mkdir(cp)
    end
    output_plot3(nt,x1,x2,4,t,size,value_x,value_y,value_z,cp,"");
end