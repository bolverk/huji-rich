function [X,Y,Pressure,Density,Vx,Vy,time,Tracers,xproc,yproc]=read_hdfMPI(filedir,filename,nproc)

filename=strcat(filedir,filename);
display('hello')
[X,Y,Pressure,Density,Vx,Vy,~,time,Tracers,~]=read_hdf(strcat(filename,sprintf('_%d.h5',0)));

xproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_x_coordinate');
yproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_y_coordinate');
NumberOfTracers=h5read(strcat(filename,sprintf('_%d.h5',0)),'/Number of tracers');
time=h5read(strcat(filename,sprintf('_%d.h5',0)),'/time');

npoints=0;
for i=0:nproc-1
   fname=strcat(filename,sprintf('_%d.h5',i));
   a=h5info(fname);
   npoints=npoints+a.Datasets(8).ChunkSize;
end
npoints
X=zeros(npoints,1);
Y=X;
Pressure=X;
Density=X;
Vx=X;
Vy=X;
if(NumberOfTracers>0)
	Tracers=zeros(npoints,NumberOfTracers);
else
	Tracers=[];
end

temp=0;
for i=0:nproc-1
    fname=strcat(filename,sprintf('_%d.h5',i));
    [Xt,Yt,Pressuret,Densityt,Vxt,Vyt,~,~,Tracerst,~]=read_hdf(fname);
    n=length(Xt);
    X(temp+1:temp+n)=Xt;
    Y(temp+1:temp+n)=Yt;
    Vx(temp+1:temp+n)=Vyt;
    Vy(temp+1:temp+n)=Vxt;
    Pressure(temp+1:temp+n)=Pressuret;
    Density(temp+1:temp+n)=Densityt;
    if(NumberOfTracers>0)
        Tracers(temp+1:temp+n,:)=Tracerst;
    end
    temp=temp+n;
end

end
