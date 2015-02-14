function [X,Y,Pressure,Density,Vx,Vy,time,Tracers,xproc,yproc]=read_hdfMPI(filedir,filename,nproc)

cd(filedir)

[X,Y,Pressure,Density,Vx,Vy,~,time,Tracers,~]=read_hdf(strcat(filename,sprintf('_%d.h5',0)));

xproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_x_coordinate');
yproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_y_coordinate');
for i=1:nproc-1
   [Xt,Yt,Pressuret,Densityt,Vxt,Vyt,~,~,Tracerst,~]=read_hdf(strcat(filename,sprintf('_%d.h5',i)));
   X=[X;Xt];
   Y=[Y;Yt];
   Pressure=[Pressure; Pressuret];
   Density=[Density;Densityt];
   Vx=[Vx;Vxt];
   Vy=[Vy;Vyt];
   Tracers=[Tracers;Tracerst];
end

end
