function [B]=func_bubblevelocity(B, tstep, D, minbubbledia_vel, ylim1, ylim2)

% ----------------------------------------------------------------------
% this function adds velocity components to the matrix B
% note: bubble velocity in frame i is 0 if
% (a) total number of bubbles in frame i is not equal to total number of bubbles in frame i+1 (coalescence/splitting/eruption)
% (b) computed vx and vy are physically unreasonable 

% Smaller bubbles distort numbering; to improve linking, use larger values for minbubbledia
% increasing minbubledia_vel and choosing [ylim1, ylim2] to exclude small bubbles may improve linking 

% max bubble velocity constraints are based on observations in visualization
% the defult setting for vxmax is that a bubble can travel atmost Width/10 in one time-step
% vymax based on bubble diameter 
% -----------------------------------------------------------------------

Btemp = B;  

% B = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR1]
nframe1 = min(B(:,1)); 
B(:,1) = B(:,1) - nframe1 + 1;                                      % modifying starting frame to begin from 1
nframe = max(B(:,1)); 
m = length(B(:,1)); B = [linspace(1,m,m)',B]; 

% B = [orig-bubble#, frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR1]
TF = B(:,8)>ylim1 & B(:,9)<ylim2 & B(:,5)>minbubbledia_vel; 
B = B(TF,:);                                                        % consider bubbles only with axial limits 
frame_nbubbles = histc(B(:,2),1:nframe);                            % # of bubbles in frame 
frame_linked = frame_nbubbles == circshift(frame_nbubbles, -1);     % equal # bubbles in consecutive frames 

B(:,6:10) = [];                             % B = [orig-bubble#, frame#, xmean, ymean, bubble-dia) 
m = length(B(:,1)); s1 = zeros(m,1); 
B = [linspace(1,m,m)', B, s1, s1];          % B = [new-bubble#, orig-bubble#, frame#, xmean, ymean, bubble-dia) 
B(:,9) = frame_linked(B(:,3));              % B(:,8) = 1 if frame linked to next frame
B = [B, s1]; 
B(:,10) = B(:,1)+frame_nbubbles(B(:,3));    % note: we don't care about non-linked frames 

% B = [new-bubble#, orig-bubble#, frame#, xmean, ymean, bubble-dia, vx, vy, linkedframe, linkedbubble] 
for i=1:length(B(:,1))   
    if B(i,9)>0 && B(i,3) ~= max(B(:,3))
        B(i,7) = (B(B(i,10),4)-B(i,4))/tstep; 
        B(i,8) = (B(B(i,10),5)-B(i,5))/tstep; 
    end
end
B(:,1) = []; 

% B = [orig-bubble#, frame#, xmean, ymean, bubble-dia, vx, vy, linkedframe, linkedbubble] 
vxmax = (D/10)/tstep; vymax = 5*0.71*sqrt(9.81.*B(:,5));
TF = abs(B(:,6))>vxmax | B(:,7)<0 | B(:,7)>vymax(:); B(TF,6:7) = 0;  
% keep only non-zero elements of B 
TF = B(:,6)==0 & B(:,7)==0; B(TF,:) = []; 

% B = [orig-bubble#, frame#, xmean, ymean, bubble-dia, vx, vy, linkedframe, linkedbubble] 

% combine matrices Btemp and B 
m=length(Btemp(:,1)); Btemp = [Btemp, zeros(m,1), zeros(m,1)]; 
for i=1:length(B(:,1))
    Btemp(B(i,1),10) = B(i,6); Btemp(B(i,1),11) = B(i,7); 
end
B = Btemp;  
% B = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]

end



