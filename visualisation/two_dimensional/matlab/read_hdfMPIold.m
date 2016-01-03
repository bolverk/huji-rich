function [X,Y,Pressure,Density,Vx,Vy,Points,time,Tracers,NumberOfPointsInCell,xproc,yproc]=read_hdfMPIold(filedir,filename,nproc,ShouldPlot,WhatToPlot,LogScale,edgestrength)
if(nargin==3),
    ShouldPlot=0;
    WhatToPlot=1;
    LogScale=0;
    edgestrength=0;
elseif (nargin==4),
    ShouldPlot=0;
    WhatToPlot=1;
    LogScale=0;
    edgestrength=0;
elseif (nargin==5),
    LogScale=0;
    edgestrength=0;
elseif (nargin==6),
    edgestrength=0;
elseif (nargin==7),
    % do nothing
else
    error('Illegal number of input arguments');
end

edgestrength=0.05;
filename=strcat(filedir,filename);
[X,Y,Pressure,Density,Vx,Vy,~,time,Tracers,~]=read_hdfold(strcat(filename,sprintf('_%d.h5',0)));

xproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_x_coordinate');
yproc=h5read(strcat(filename,sprintf('_%d.h5',0)),'/proc_y_coordinate');
h=h5info(strcat(filename,sprintf('_%d.h5',0)));
NumberOfTracers=h5read(strcat(filename,sprintf('_%d.h5',0)),'/Number of tracers');
time=h5read(strcat(filename,sprintf('_%d.h5',0)),'/time');

npoints=0;
for i=0:nproc-1
    fname=strcat(filename,sprintf('_%d.h5',i));
    a=h5info(fname);
    npoints=npoints+a.Datasets(10).ChunkSize;
end
maxfaces=0;
for i=0:nproc-1
    fname=strcat(filename,sprintf('_%d.h5',i));
    nVert=h5read(fname,'/Number of vertices in cell');
    maxfaces=max(max(nVert),maxfaces);
    if(maxfaces>20)
        display('Warning, max number of faces exceeds 20!!')
    end
end
npoints
X=zeros(npoints,1);
Y=X;
Pressure=X;
Density=X;
Vx=X;
Vy=X;
Points=zeros(npoints,maxfaces,2);
NumberOfPointsInCell=X;
if(NumberOfTracers>0)
    Tracers=zeros(npoints,NumberOfTracers);
else
    Tracers=[];
end
temp=0;
for i=0:nproc-1
    fname=strcat(filename,sprintf('_%d.h5',i));
    [Xt,Yt,Pressuret,Densityt,Vxt,Vyt,Pointst,~,Tracerst,nincellt]=read_hdfold(fname);
    n=length(Xt);
    X(temp+1:temp+n)=Xt;
    Y(temp+1:temp+n)=Yt;
    Vx(temp+1:temp+n)=Vxt;
    Vy(temp+1:temp+n)=Vyt;
    Pressure(temp+1:temp+n)=Pressuret;
    Density(temp+1:temp+n)=Densityt;
    temp3=size(Pointst);
    temp2=zeros(temp3(1),maxfaces,2);
    temp2(:,1:temp3(2),:)=Pointst;
    Points(temp+1:temp+n,:,:)=temp2;
    NumberOfPointsInCell(temp+1:temp+n)=nincellt;
    if(NumberOfTracers>0)
        Tracers(temp+1:temp+n,:)=Tracerst;
    end
    temp=temp+n;
end
display('Finished reading data');
if(ShouldPlot==1||ShouldPlot==2)
    NumberOfCells=length(X);
    Vertices=zeros(maxfaces*NumberOfCells,2);
    Faces=zeros(maxfaces,NumberOfCells);
    for i=1:NumberOfCells
        n=NumberOfPointsInCell(i);
        if(ShouldPlot==1||ShouldPlot==2)
            Vertices((i-1)*maxfaces+1:(i-1)*maxfaces+n,1)=Points(i,1:n,1);
            Vertices((i-1)*maxfaces+1:(i-1)*maxfaces+n,2)=Points(i,1:n,2);
            Faces(1:n,i)=((i-1)*maxfaces+1):((i-1)*maxfaces+n);
            Faces(n+1:maxfaces,i)=NaN;
        end
    end
end

if(ShouldPlot==1||ShouldPlot==2)
    if(ShouldPlot==1)
        f1=figure;
        set(f1,'Units','normalized')
        set(f1, 'Position', [0.03 0.03 0.65 0.85])
    end
    hold on;
    maxdraw=7500;
    nloops=ceil(length(Density)/maxdraw);
    for i=1:nloops
        maxindex=min(i*maxdraw,length(Density));
        switch (WhatToPlot)
            case 1
                if(LogScale==1)
                    caxis([min(log10(Density)) max(log10(Density*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',log10(Density((i-1)*maxdraw+1:maxindex)),'FaceColor','flat','EdgeAlpha',edgestrength);
                else
                    caxis([min(Density) max(Density)*1.01]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',Density((i-1)*maxdraw+1:maxindex),'FaceColor','flat','EdgeAlpha',edgestrength);
                end
            case 2
                if(LogScale==1)
                    caxis([min(log10(Pressure)) max(log10(Pressure*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',log10(Pressure((i-1)*maxdraw+1:maxindex)),'FaceColor','flat','EdgeAlpha',edgestrength);
                else
                    caxis([min(Pressure) max(Pressure)*1.01]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',Pressure((i-1)*maxdraw+1:maxindex),'FaceColor','flat','EdgeAlpha',edgestrength);
                end
            case 3
                if(LogScale==1)
                    caxis([min(log10(Temperature)) max(log10(Temperature*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',log10(Temperature((i-1)*maxdraw+1:maxindex)),'FaceColor','flat','EdgeAlpha',edgestrength);
                else
                    caxis([min(Temperature) max(Temperature)*1.01]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',Temperature((i-1)*maxdraw+1:maxindex),'FaceColor','flat','EdgeAlpha',edgestrength);
                end
            case 4
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Change here gamma as needed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                gamma=5/3;
                Entropy=Pressure./Density.^gamma;
                if(LogScale==1)
                    caxis([min(log10(Entropy)) max(log10(Entropy*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',log10(Entropy((i-1)*maxdraw+1:maxindex)),'FaceColor','flat','EdgeAlpha',0);
                else
                    caxis([min(Entropy) max(Entropy)*1.01]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',Entropy((i-1)*maxdraw+1:maxindex),'FaceColor','flat','EdgeAlpha',0);
                end
            case 5
                tracerindex=1;
                if(LogScale==1)
                    caxis([min(log10(Tracers(:,tracerindex))) max(log10(Tracers(:,tracerindex)*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',log10(Tracers((i-1)*maxdraw+1:maxindex,tracerindex)'),'FaceColor','flat','EdgeAlpha',0);
                else
                    caxis([min((Tracers(:,tracerindex))) max((Tracers(:,tracerindex)*1.01))]);
                    patch('Faces',Faces(:,(i-1)*maxdraw+1:maxindex)','Vertices',Vertices,'FaceVertexCData',Tracers((i-1)*maxdraw+1:maxindex,tracerindex)','FaceColor','flat','EdgeAlpha',0.05);
                end
        end
    end
    colorbar;
    axis equal;
end

end
