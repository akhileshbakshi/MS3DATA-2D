function [nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2)
 
bubblepropertiestotal = [0 0 0 0 0 0 0 0 0];

fileID = fopen(bubblefile);
Anet = textscan(fileID, '%f %f %f');
fclose(fileID);
Anet = [Anet{1,1} Anet{1,2} Anet{1,3}];
 
% add fictitious row in A 
Anet = [Anet; 10000, 10000, 10000];
[s1,s2] = ismember(unique(Anet(:,1)),Anet(:,1));
if nframes==0
    nframes = length(unique(Anet(:,1)))-1;                  % to account for the extra 10000 at the end 
end
frameloc = s2; 

[R,nx,H,ny,coarsegridglobal] = func_readgeometry;    % function call is expensive in parfor
% coarsegridglobal = [cellgeom, xgeom, ygeom];

parfor framei = 1:nframes

framei

% --------------------------------------------------------------
% If single file used, use this snippet 
Anet_local = Anet;
frameloc_local = frameloc; 

nz =1;                                      % For Cartesian 2D 
A1 = Anet(frameloc(framei):frameloc(framei+1)-1,2);     
A2 = Anet(frameloc(framei):frameloc(framei+1)-1,3);  
epgcoarse = epgcutoff*ones(nx*ny*nz,1); 
epgcoarse(A1) = A2;                         % replacing epg where epg>epgcutoff

% coarsegrid = [coarsegrid epgcoarse];      % coarsegrid = [cell#, x, y, epg]
coarsegrid = coarsegridglobal; 
xcoarse = coarsegrid(:,2); 
ycoarse = coarsegrid(:,3); 

xgrid = linspace(0,R,xsmooth*nx); 
ygrid = linspace(0,H,ysmooth*ny);
[xgrid,ygrid] = meshgrid(xgrid,ygrid);
xgrid = reshape(xgrid,[],1);
ygrid = reshape(ygrid,[],1);
 
epggrid = griddata(xcoarse,ycoarse,epgcoarse,xgrid,ygrid);


% ----------------------------------------------------------------
% restructuring matrix

deltax = R/(xsmooth*nx-1);          % in fine grid 
deltay = H/(ysmooth*ny-1); 

B = [xgrid ,ygrid, epggrid]; 
B(isnan(B)) = 0;        % boundary cells have NaN (no interpolation) 
B = sortrows(B,2);      % sortrpws based on axial location 

% clear xgrid; clear ygrid; clear epggrid;

% ----------------------------------------------------------------
% checking neighbouring cells for linking 
% B = [# x y epg]

m = length(B(:,1)); 
v = zeros(m,1);          
B = [linspace(1,m,m)',B, v, v];   % numbering & 2 columns to store neighbours

% B = [# x y epg neighbor1 neighbor2]
B(:,5) = B(:,4)>epgbubble & circshift(B(:,4),-1)>epgbubble;                % right neighbour
B(:,5) = (B(:,1)+1).* B(:,5);

B(:,6) = B(:,4)>epgbubble & circshift(B(:,4),-xsmooth*nx)>epgbubble;       % top neighbour
B(:,6) = (B(:,1)+xsmooth*nx).* B(:,6);

TF = B(:,5)==0; 
s = TF.*B(:,6); 
B(:,5)=B(:,5)+s; 
B(:,6)=B(:,6)-s;


% ----------------------------------------------------------------
% checking if bubble is at interface (1=interface 0=interior)
% B = [# x y z epg neighbour#1 neighbour#2 neighbour#3]

v = zeros(length(B(:,1)),1); 
B = [B, v]; 

s = B(:,4);
TF1 = s>epgbubble;
TF2 = circshift(s,1)>epgbubble;
TF3 = circshift(s,-1)>epgbubble;
TF4 = circshift(s,xsmooth*nx)>epgbubble;
TF5 = circshift(s,-xsmooth*nx)>epgbubble;
B(:,7) = TF1 & TF2 & TF3 & TF4 & TF5; 
B(:,7) = 1-B(:,7); 

TF = B(:,4)>epgbubble; 
B = B(TF,:); 

% ----------------------------------------------------------------
% renumbering indices of points 
% B = [# x y epg neighbour#1 neighbour#2 interface]

[s1,s2] = ismember(B(:,5),B(:,1)); B(:,5) = s2; 
[s1,s2] = ismember(B(:,6),B(:,1)); B(:,6) = s2; 

B(:,1) = []; 
B(:,3) = [];   


% ----------------------------------------------------------------
% linking bubbles
% B = [x y neighbour#1 neighbour#2 interface]

v = zeros(length(B(:,1)),1); B = [B, v, v];           
s = B(:,5); B(:,5)=[]; B = [B, s];
% B = [x y neighbour#1 neighbour#2 bubble#1 bubble#2 interface]                         

% start numbering bubbles 
bubblenumctr = 0;           % keeps note of latest bubblenum

for i=1:length(B(:,1)) 
    if B(i,5)>0 
        if B(i,3)>0 
            if B(B(i,3),5) ==0 
                B(B(i,3),5)= B(i,5); 
            elseif B(B(i,3),5)~=B(i,5)
                B(B(i,3),6)=B(i,5); end 
            if B(i,4)>0                     % because B(i,4)>0 only possible if B(i,3)>0
                B(B(i,4),5)= B(i,5); end 
        end
        
    elseif B(i,3)>0 && B(B(i,3),5)>0        % B(i,5) = 0, and B(i,4) will not be numbered
            B(i,5) = B(B(i,3),5); 
        if B(i,4)>0
            B(B(i,4),5)= B(i,5); end     
        
    else                                    % B(i,5)=0 and (B(i,3)=0 || B(B(i,3),5)=0)
        bubblenumctr = bubblenumctr+1; 
        B(i,5) = bubblenumctr; 
        if B(i,3)>0      
            B(B(i,3),5)= B(i,5); end
        if B(i,4)>0
            B(B(i,4),5)= B(i,5); end     
    end 
end 


% ----------------------------------------------------------------
% find unique dispute combinations (i.e. B(i,5)~=B(i,6))
% B = [x y neighbour#1 neighbour#2 bubble#1 bubble#2 interface]                         

dispute = [B(:,5:6), B(:,6) > 0]; 
% dispute = [bubble#1 bubble#2 bubble#3 dispute]
TF = dispute(:,end)<1; dispute(TF,:)=[];    % keep only dispute cases

% get unqiue ID of combination (assume max 9999 disputes allowed) 

dispute(:,3) = []; 
dispute = sort(dispute,2,'descend'); 
dispute(:,3) = dispute(:,1)*10^4 + dispute(:,2);
disputecombination = zeros(length(unique(dispute(:,3))),3); 
disputecombination(:,1) = mod(unique(dispute(:,3)),10^4);
disputecombination(:,2) = (unique(dispute(:,3))-disputecombination(:,1))/10^4; 
disputecombination=sortrows(disputecombination,2); 

% clear dispute; 
% checking for possibility that c->a and c->b 

[m n] = size(disputecombination);

if m>1          % consider cases with >1 disputes only 

    % add a fictitious row in dipute combination for cases where only two entries [a b; a c] 
    disputecombination = [zeros(1,3); disputecombination];
    s = disputecombination(:,2); 
    % disputecombination = [disputecombination, zeros(length(s),1)];
    disputecombination(:,3) = s == circshift(s,1);
    conflicts = sum(disputecombination(:,3));
    disputecombination(1,:) = []; 

    while conflicts>0  
        [s1, s2] = ismember(1,disputecombination(:,3));
        disputecombination(s2,2)= max(disputecombination(s2,1),disputecombination(s2-1,1));
        disputecombination(s2,1)= min(disputecombination(s2,1),disputecombination(s2-1,1));

        s1 = disputecombination(:,1)*10^4 + disputecombination(:,2);

        disputecombinationtemp = zeros(length(disputecombination(:,1)),3); 
        disputecombinationtemp(1:length(unique(s1)),2) = mod(unique(s1),10^4);
        disputecombinationtemp(1:length(unique(s1)),1) = (unique(s1)-disputecombinationtemp(1:length(unique(s1)),2))/10^4; 
        if (disputecombinationtemp(end,1)==0) 
            disputecombinationtemp(end,:)=[];
        end           
        disputecombinationtemp=sortrows(disputecombinationtemp,2); 

        s = disputecombinationtemp(:,2); 
        disputecombinationtemp(:,3) = s == circshift(s,1);
        conflicts = sum(disputecombinationtemp(:,3));
        disputecombination = disputecombinationtemp; 
        disputecombinationtemp = zeros(length(disputecombination(:,1)),3);
    end
end 

% ----------------------------------------------------------------
% modify bubble# for disputed cases
% B = [x y neighbour#1 neighbour#2 bubble#1 bubble#2 interface]
% dispute combination = [newvalue oldvalue conflict] 


% change B from highest dispute case e.g. 19->10 & 10->5 => 19->5 
 
for i=length(disputecombination(:,1)):-1:1
    s1 = B(:,5);
    s1(s1==disputecombination(i,2))=disputecombination(i,1);
    numbermod = find(s1>disputecombination(i,2));            % renumbering the higher bubble    
    s1(numbermod) = s1(numbermod)-1; 
    B(:,5) = s1;  
end 

B(:,6)=[];
B = sortrows(B,5);
B(:,3:4) = [];   

% B = [x y bubble#1 interface]   

bubbleproperties = func_bubbleproperties(B(:,1:3), deltax, deltay); 
% bubbleproperties = [bubble#, xmean, ymean, Area, xmin, xmax, ymin, ymax, AR]
bubbleproperties(:,4) = (4*bubbleproperties(:,4)/pi).^0.5; 
TF1 = bubbleproperties(:,8)-bubbleproperties(:,7)<mincordlength; 
TF2 = bubbleproperties(:,6)-bubbleproperties(:,5)<minCSlength;
TF3 = bubbleproperties(:,4)<minbubbledia; 
bubbleproperties(TF1 | TF2 | TF3,:) = []; 

% the first column now gets frame number 
bubbleproperties(:,1) = framei;
bubblepropertiestotal = [bubblepropertiestotal; bubbleproperties];
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR]

end
bubblepropertiestotal(1,:) = []; 

TF1 = bubblepropertiestotal(:,8) > ycutoff2 | bubblepropertiestotal(:,7) < ycutoff1 ; 
bubblepropertiestotal(TF1,:) = [];        % remove bubbles touching top and bottom of frame      

end







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






function [C] = func_bubbleproperties(B, deltax, deltay)

% B = [x y bubble#]

% Adding fictitious row in B for ease of indices later 
B = [B; 10000, 10000, 10000]; 
[s1,s2] = ismember(unique(B(:,3)),B(:,3));
C = [unique(B(:,3)), s2];           % bubble number and first cell 

% last row in C will be fictitious (due to fictitious last row in B) 
for i=1:size(C(:,1))-1
    C(i,3) = mean(B(C(i,2):C(i+1,2)-1,1));      % xmean 
    C(i,4) = mean(B(C(i,2):C(i+1,2)-1,2));      % ymean
    C(i,5) = (C(i+1,2)-C(i,2))*deltax*deltay;   % Area       
    C(i,6) = min(B(C(i,2):C(i+1,2)-1,1));       % xmin
    C(i,7) = max(B(C(i,2):C(i+1,2)-1,1));       % xmax
    C(i,8) = min(B(C(i,2):C(i+1,2)-1,2));       % ymin
    C(i,9) = max(B(C(i,2):C(i+1,2)-1,2));       % ymax
    C(i,10) =(C(i,9)-C(i,8))/(C(i,7)-C(i,6));   % AR
end

B(end,:) = []; C(end,:) = []; 
% sortrows based on ylocation 
C = sortrows(C,4); C(:,2) = []; 
% renumber bubbles 
C(:,1) = linspace(1,length(C(:,1)),length(C(:,1)))';
% output = [bubble#, xmean, ymean, bubbleA, xmin, xmax, ymin, ymax, AR]

end


