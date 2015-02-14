function [X,Pressure,Density,Vx,time]=read_hdf1D(filename,ShouldPlot,WhatToPlot,LogScale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab script to read RICH binary files in float/double format
%
%   Inputs  : filename - The location of the binary file 
%             ShouldPlot - 1 indicates that a plot is wanted 0 otherwise
%             WhatToPlot - 1: density, 2: pressure, 3: pressure/density, 4:
%             Entropy (assuming gamma=5/3)
%             LogScale - 1 Plot Log scaled, 0 linear scale
%   Output  : X - The mesh points
%             Pressure - The pressure
%             Density - The density
%             Vx - The x component of the velocity
%             time - The time of the simulation
%   Example : [X,Pressure,Density,xVelocity,time]=read_hdf1D("output.bin",1,1,1)
%               Plots a log scaled density plot from the file "output.bin"
%             [X,Pressure,Density,xVelocity,time]=read_hdf1D("output.bin",0)
%               Only reads the data from "output.bin"
%             [X,Pressure,Density,xVelocity,time]=read_hdf1D("output.bin")
%               Only reads the data from "output.bin"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin==1),
ShouldPlot=0;
WhatToPlot=1;
LogScale=0;
elseif (nargin==2),
  ShouldPlot=0;
  WhatToPlot=1;
  LogScale=0;
elseif (nargin==3),
   LogScale=0;
elseif (nargin==4),
   % do nothing
else
   error('Illigal number of input arguments');
end

% Read the HDF5 data
Density=h5read(filename,'/density');
Pressure=h5read(filename,'/pressure');
X=h5read(filename,'/grid');
Vx=h5read(filename,'/x_velocity');
time=h5read(filename,'/time');


draw=WhatToPlot;
Log=LogScale;
Temperature=Pressure./Density;


if(ShouldPlot==1)
    f1=figure;
    set(f1,'Units','normalized')
    set(f1, 'Position', [0.03 0.03 0.65 0.85])
    hold on;
    switch (draw)
        case 1
            if(Log==1)
                plot(X,log10(Density),'x');
            else
                plot(X,Density,'x');
            end
        case 2
            if(Log==1)
                plot(X,log10(Pressure),'x');
            else
                plot(X,Pressure,'x');            
            end
        case 3
            if(Log==1)
                plot(X,log10(Temperature),'x');
            else
                plot(X,Temperature,'x');
            end
        case 4
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Change here gamma as needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gamma=5/3;
            Entropy=Pressure./Density.^gamma;
            if(Log==1)
                plot(X,log10(Entropy),'x');
            else
                plot(X,Entropy,'x');
            end
    end
end
