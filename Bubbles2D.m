clear all; clc; 

% 1. Mfix file properties 
nframes = 0;            % 0 to read from void fraction data file
                        % non-zero to specify exact number of frames to perform stats
epgcutoff = 0.65; 
ycutoff2   = 0.69;      % maximum domain y-extremity (ycutoff2 must be < ymax from simulation data) 
ycutoff1 =  0.004;      % minimum domain y-extremity
                        % domain extremeties are required to remove ambiguous bubbles (touching freeboard and distributor) 

% 2. Modify Geometry.xlsx and enter other simulation data
D = 0.3;                % this is the bed diameter (kept R to preserve nomenclature with cylindrical 3D) 
tstep = 0.01;           % indicates frequency of void fraction data 

% 3. input/output files names 
% sample file provided has data corresponding to 300 frames 
bubblefile = 'bubblestats2D.txt';
printfile = 'bubbles2D'; 

% 4. bubble properties  
epgbubble = 0.7;        % identifies bubbles
mincordlength = 0.01;   % discard bubbles which are very small  
minCSlength = 0.01;     % to avoid extremely small bubbles (infinite AR)
minbubbledia = 0.01;    % discard bubbles which are very small  
                        % NOTE: discarding small bubbles helps bubbling linking 

% 5. Grid smoothening 
ysmooth = 1;            % grid is refined based on this factor 
xsmooth = 1;           

% 5. limits for post processing for isolating bubbles within region of interest
ylim1 = 0;  
ylim2 = ycutoff2; 
rlim1 = 0; 
rlim2 = D;             

% 6. Statistics for average computations 
nbinsax = 15;
nbinsrad = 4; 
nbinspdf = 20; 

% ----------------------------------------------------------------
[nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2);
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR]

bubblepropertiestotal = func_bubblevelocity(bubblepropertiestotal, tstep, D); 
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]
% note: bubble velocity in frame i is 0 if
% (a) total number of bubbles in frame i is not equal to total number of bubbles in frame i+1 (coalescence/splitting/eruption)
% (b) computed vx and vy are physically unreasonable 

% ----------------------------------------------------------------
% sample for computing average statistics 
[bubblestats_2D, bubblestats_ax, bubblestats_rad]=func_bubblestatistics(bubblepropertiestotal, nframes, nbinsax, nbinsrad, ylim1,ylim2,rlim1,rlim2);
% bubblestats_2D = [binr, biny, nb, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 
% bubblestats_ax = [biny, nb_y, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 
% bubblestats_rad= [binr, nb_r, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 

% % ----------------------------------------------------------------
% % sample for writing to files 
% filename = strcat(printfile,'_BubbleStats_Ax.txt');
% dlmwrite(filename,bubblestats_ax,'delimiter',' ','precision',4); 







