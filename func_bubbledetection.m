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

