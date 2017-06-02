% variables are tailored for 'bubblestats2D.txt' and 'geometry.xlsx' provided 
% Major Edit 1: func_bubblevelocity has been edited to include minbubbledia_vel, ylim1, ylim2 for better linking 
% Major Edit 2: option lagrangetracking has been added where bubbles are individually tracked through frames 

clear all; clc; 

% 1. Mfix file properties 
nframes = 0;            % 0 to read from void fraction data file
                        % non-zero to specify exact number of frames 
epgcutoff = 0.65;       % consider only cells so that ep_g>epgcutoff; epgcutoff<epgbubble
ycutoff2   = 0.69;      % maximum domain y-extremity (ycutoff2 must be < ymax from simulation data) 
ycutoff1 =  0.004;      % minimum domain y-extremity
                        % domain extremeties are required to remove ambiguous bubbles (touching freeboard and distributor) 

% 2. Modify Geometry.xlsx and enter other simulation data
D = 0.3;                % diameter/width of bed/slice
tstep = 0.01;           % time step of data sampling 

% 3. input/output files names 
% sample file provided has data corresponding to 300 frames 
bubblefile = 'bubblestats2D.txt';
printfile = 'bubbles2D'; 

% 4. criteria for bubble detection   
epgbubble = 0.7;        % threshold voidage for bubble (interphase) detection 
mincordlength = 0.01;   % discard small bubbles   
minCSlength = 0.01;     % discard small bubbles 
minbubbledia = 0.01;    % discard small bubbles  
ysmooth = 1;            % y-grid refinement
xsmooth = 1;            % x-grid refinement          

% 5. criteria for postprocessing of detected bubbles
ylim1 = 0;              % min y for postprocessing 
ylim2 = ycutoff2;       % max y for postprocessing
rlim1 = 0;              % min x for postprocessing
rlim2 = D;              % max x for postprocessing
minbubbledia_vel = 0.01;% discard very small bubbles for bubble linking
diaratio = 1.1;         % maximum permissible ratio of bubble dia for linking  
dmax = 0.05;            % maximum permissible distance traveled by bubble in one time-step  
tolerance  = 0.0;       % minimum permissible bubble y-velocity = -tolerance x time-step 
lagrangetracking = 1;   % (recommended) set 1 to turn on lagrangian tracking of bubbles 
                        % linking is affected by bubble activity- splitting, coalescence and eruption
                        % for best linking results write data at high frequency 
                        % if globaltracking (lagrangetracking=0), consider increasing minbubledia_vel
                        % and choosing [ylim1, ylim2] to exclude small bubbles may improve linking 

% 6. Statistics for average computations 
nbinsax = 10;           % # bins for axial statistics between [ylim1, ylim2]
nbinsrad = 4;           % # bins for radial/lateral statistics between [rlim1, rlim2]

% ----------------------------------------------------------------
[nframes, bubblepropertiestotal] = func_bubbledetection(bubblefile, xsmooth, ysmooth, epgcutoff, epgbubble, mincordlength, minCSlength, minbubbledia, nframes, ycutoff1, ycutoff2);
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR]

[bubbletrace, bubblepropertiestotal] = func_bubblevelocity(bubblepropertiestotal, tstep, minbubbledia_vel, ylim1, ylim2, lagrangetracking, diaratio, dmax, tolerance); 
% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]

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

% ----------------------------------------------------------------
% sample for plotting bubbletrace- only possible if lagrangetracking = 1

tracebubblenum = [9 22 161 165];
for j=1:length(tracebubblenum)
  switch j 
    case 1; color = 'r';
    case 2; color = 'k'; 
    case 3; color = 'b';
    case 4; color = 'm';
  end 

  traceindex = find(bubbletrace(:,tracebubblenum(j))); 
  tracebubbles = bubbletrace(traceindex,tracebubblenum(j)); 
  xscatter = bubblepropertiestotal(tracebubbles,2); yscatter = bubblepropertiestotal(tracebubbles,3); size = bubblepropertiestotal(tracebubbles,4).^2; 
  scatter(xscatter,yscatter,500000*size,color,'LineWidth',1.2);
  xlabel ('x [m]','FontWeight','bold','fontsize',20);
  ylabel ('y [m]','FontWeight','bold','fontsize',20);
  xlim([0 0.3]); ylim([0 0.6]);
  hold on;
end








