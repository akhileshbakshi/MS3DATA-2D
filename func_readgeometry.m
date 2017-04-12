function [R,nr,H,ny, A] = func_readgeometry()

% Note: This function uses for loops which are potentially slow but since
% this function only needs to be called once, it is okay ! 

File = 'Geometry.xlsx';
R = xlsread(File,'Sheet1','C3');
nr = xlsread(File,'Sheet1','C4');
H = xlsread(File,'Sheet1','E3');
ny = xlsread(File,'Sheet1','E4'); 
nz = 1; 

% x direction 
[num,txt,raw] = xlsread(File,'Sheet1','C5');
if strcmp(txt,'no')
    drrange = strcat('C7:C',num2str(7+nr-1));
    dr = xlsread(File,'Sheet1',drrange);
else dr = (R/nr)*ones(nr,1);
end

% y direction 
[num,txt,raw] = xlsread(File,'Sheet1','E5');
if strcmp(txt,'no')
    dyrange = strcat('C7:C',num2str(7+ny-1));
    dy = xlsread(File,'Sheet1',drrange);
else dy = (H/ny)*ones(ny,1);
end

% setup coarse mesh in x direction 
xtemp(1)=dr(1);         % to get wall centers 
for i=2:size(dr)
    xtemp(i)=0;
    for k=1:i
        xtemp(i)=xtemp(i)+dr(k); end
end
xtemp = [0 xtemp];
for i=1:size(dr)        % to get cell centers 
    xcell(i)=0.5*(xtemp(i+1)+xtemp(i)); end

% setup coarse grid in y direction 
ytemp(1)=dy(1);         % to get wall centers 
for i=2:size(dy)
    ytemp(i)=0;
    for k=1:i
        ytemp(i)=ytemp(i)+dy(k);
    end
end
ytemp = [0 ytemp];
for i=1:size(dy)        % to get cell centers 
    ycell(i)=0.5*(ytemp(i+1)+ytemp(i));
end

xgeom = zeros(nr*ny*nz,1);
ygeom = zeros(nr*ny*nz,1);
cellgeom = linspace(1,nr*ny*nz,nr*ny*nz); 

% the cell ordering must replicate IJK in post-mfix
% assumes r first, then y, and then theta 

for i=1:nr*ny 
    if mod(i,nr)>0
        xgeom(i)= xcell(mod(i,nr));
        ygeom(i)= ycell((i-mod(i,nr))/nr+1);
    else
        xgeom(i) = xcell(end);
        ygeom(i) = ycell(i/nr);
    end   
end

cellgeom = cellgeom'; 

A = [cellgeom, xgeom, ygeom];

clear xtemp; clear ytemp; clear xcell; clear ycell; clear cellgeom; 

end

