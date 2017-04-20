% variables are tailored for 'bubblestats2D.txt' and 'geometry.xlsx' provided 
% Major Edit 1: func_bubblevelocity has been edited to include minbubbledia, ylim1, ylim2 for better linking 

clear all; clc; 

% 1. Mfix file properties 
nframes = 0;            % 0 to read from void fraction data file
                        % non-zero to specify exact number of frames to perform stats
epgcutoff = 0.65; 
ycutoff2   = 0.69;      % maximum domain y-extremity (ycutoff2 must be < ymax from simulation data) 
ycutoff1 =  0.004;      % minimum domain y-extremity
                        % domain extremeties are required to remove ambiguous bubbles (touching freeboard and distributor) 

% 2. Modify Geometry.xlsx and enter other simulation data
D = 0.3;                % bed diameter 
tstep = 0.01;           % indicates frequency of void fraction data 

% 3. input/output files names 
% sample file provided has data corresponding to 300 frames 
bubblefile = 'bubblestats2D.txt';
printfile = 'bubbles2D'; 

% 4. criteria for bubble detection   
epgbubble = 0.7;        % identifies bubbles
mincordlength = 0.01;   % discard bubbles which are very small  
minCSlength = 0.01;     % to avoid extremely small bubbles (infinite AR)
minbubbledia = 0.01;    % discard bubbles which are very small  
                        % NOTE: discarding small bubbles helps bubbling linking 

% 5. Grid smoothening 
ysmooth = 1;            % grid is refined based on this factor 
xsmooth = 1;           

% 5. criteria for postprocessing of detected bubbles
ylim1 = 0;  
ylim2 = ycutoff2; 
rlim1 = 0; 
rlim2 = D;  
minbubbledia_vel = 0.02; 

% 6. Statistics for average computations 
nbinsax = 10;
nbinsrad = 4; 

% ----------------------------------------------------------------
[nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2);
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR]

bubblepropertiestotal = func_bubblevelocity(bubblepropertiestotal, tstep, D, minbubbledia_vel, ylim1, ylim2); 
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]
% note: bubble velocity in frame i is 0 if
% (a) total number of bubbles in frame i is not equal to total number of bubbles in frame i+1 (coalescence/splitting/eruption)
% (b) computed vx and vy are physically unreasonable 
% increasing minbubledia_vel and choosing [ylim1, ylim2] to exclude small bubbles may improve linking 

% ----------------------------------------------------------------
% sample for computing average statistics 
[bubblestats_2D, bubblestats_ax, bubblestats_rad]=func_bubblestatistics(bubblepropertiestotal, nbinsax, nbinsrad, ylim1,ylim2,rlim1,rlim2);
% bubblestats_2D = [binr, biny, nb, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 
% bubblestats_ax = [biny, nb_y, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 
% bubblestats_rad= [binr, nb_r, area-dia, CSmax, cord, AR, nbubbles_linked, abs(vx), vy]; 

% % ----------------------------------------------------------------
% % sample for writing to files 
% filename = strcat(printfile,'_BubbleStats_Ax.txt');
% dlmwrite(filename,bubblestats_ax,'delimiter',' ','precision',4); 







