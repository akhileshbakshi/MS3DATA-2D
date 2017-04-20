function [B]=func_bubblevelocity(B, tstep, D)

Btemp = B;  

% B = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR1]

nframe1 = min(B(:,1)); 
B(:,1) = B(:,1) - nframe1 + 1;                                      % modifying starting frame to begin from 1
frame_nbubbles = histc(B(:,1),1:max(B(:,1)));                       % # of bubbles in frame 
frame_linked = frame_nbubbles == circshift(frame_nbubbles, -1);     % equal # bubbles in consecutive frames 

B(:,5:9) = [];                              % B = [frame#, xmean, ymean, area) 
s1 = zeros(length(B(:,1)),1); 
B = [B, s1, s1]; 
B(:,7) = frame_linked(B(:,1));              % B(:,7) = 1 if frame linked to next frame
B = [linspace(1,length(B(:,1)),length(B(:,1)))',B, s1]; 
% B = [bubble#, frame#, xmean, ymean, bubble-dia, vx, vy, linkedframe, linkedbubble] 

B(:,9) = B(:,1)+frame_nbubbles(B(:,2));     % note: we don't care about non-linked frames 

for i=1:length(B(:,1))   
    if B(i,8)>0 && B(i,2) ~= max(B(:,2))
        B(i,6) = (B(B(i,9),3)-B(i,3))/tstep; 
        B(i,7) = (B(B(i,9),4)-B(i,4))/tstep; 
    end
end

% B = [bubble#, frame#, xmean, ymean, bubble-dia, vx, vy, linkedframe, linkedbubble] 
vxmax = (D/10)/tstep;      
vymax = 3*0.71*sqrt(9.81.*B(:,5));
TF = abs(B(:,6))>vxmax | B(:,7)<0 | B(:,7)>vymax(:);  
B(TF,6:7) = 0;  

% since no row has been deleted, number of rows and order (bubbles) in Btemp and B are preserved 
B = [Btemp B(:,6:7)]; 
% B = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]

end



