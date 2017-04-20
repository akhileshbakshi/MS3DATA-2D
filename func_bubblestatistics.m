function [bubblestats_2D, bubblestats_ax, bubblestats_rad]=func_bubblestatistics(B, nbinsax, nbinsrad, y1,y2, r1, r2)

% bubblepropertiestotal = [frame#, xmean, ymean, bubble-dia, xmin, xmax, ymin, ymax, AR, vx, vy]
B(:,4) = 0.25*pi*B(:,4).^2; 

% bubblepropertiestotal = [frame#, xmean, ymean, area, xmin, xmax, ymin, ymax, AR, vx, vy]

% remove bubbles outside the range specified 
TF = B(:,3)<y1 | B(:,3)>y2 | abs(B(:,2))<r1 |  abs(B(:,2))>r2 ; 
B(TF,:) = [];

bins_ax=linspace(y1,y2,nbinsax+1)';         % axial bin edges
bincenter_ax = 0.5*(bins_ax+circshift(bins_ax,-1));   
bincenter_ax(end) = []; 

bins_rad=linspace(r1,r2,nbinsrad+1)';       % radial bin edges
bincenter_rad = 0.5*(bins_rad+circshift(bins_rad,-1));   
bincenter_rad(end) = []; 
 
nbubbles = zeros(nbinsax*nbinsrad,1); 
bubblevolavg = zeros(nbinsax*nbinsrad,1); 
bubbleAR1avg = zeros(nbinsax*nbinsrad,1); 
bincenter_2Dr = zeros(nbinsax*nbinsrad,1); 
bincenter_2Dy = zeros(nbinsax*nbinsrad,1); 
bubbleCSavg = zeros(nbinsax*nbinsrad,1);
bubblecordavg = zeros(nbinsax*nbinsrad,1);
bubblevxavg = zeros(nbinsax*nbinsrad,1);
bubblevyavg = zeros(nbinsax*nbinsrad,1);
nbubbles_vel = zeros(nbinsax*nbinsrad,1); 

% gather bubble stats in every bin- row first and then along height
for i=1:nbinsax*nbinsrad 
    i_y = 1+(i-1-mod(i-1,nbinsrad))/nbinsrad;
    i_r = i-(i_y-1)*nbinsrad; 
      
    TF1 = B(:,3)>=bins_ax(i_y) & B(:,3)<bins_ax(i_y+1);             % check axial location
    distance_r = abs(B(:,2)); 
    TF2 = distance_r>=bins_rad(i_r) & distance_r<bins_rad(i_r+1);   % checking radial location
    TF = TF1 & TF2; 
    
    bincenter_2Dr(i) = bincenter_rad(i_r); 
    bincenter_2Dy(i) = bincenter_ax(i_y); 
    nbubbles(i) = sum(TF); 
    if nbubbles(i)>0 
        bubblevolavg(i) = mean(B(TF,4));         
        bubbleAR1avg(i) = mean(B(TF,9));
        bubbleCSavg(i) = mean(B(TF,6)-B(TF,5));
        bubblecordavg(i) = mean(B(TF,8)-B(TF,7));
    end
   
    TF3 = B(:,11)>0;                                                % ensure averaging if bubble is linked 
    TF = TF & TF3; 
    nbubbles_vel(i) = sum(TF); 
    if nbubbles_vel(i)>0 
        bubblevxavg(i) = sum(B(TF,4).*abs(B(TF,10)))/sum(B(TF,4));  % area weighted avg         
        bubblevyavg(i) = sum(B(TF,4).*abs(B(TF,11)))/sum(B(TF,4));  % area weighted avg  
    end   
end 

% aggregate radial bins for axial variation 
nbubbles_ax = zeros(nbinsax,1); 
bubblevolavg_ax = zeros(nbinsax,1); 
bubbleAR1avg_ax = zeros(nbinsax,1); 
bubbleCSavg_ax = zeros(nbinsax,1); 
bubblecordavg_ax = zeros(nbinsax,1);
bubblevxavg_ax = zeros(nbinsax,1);
bubblevyavg_ax = zeros(nbinsax,1); 
nbubbles_vel_ax = zeros(nbinsax,1); 
for i=1:nbinsax 
     i1 = 1+(i-1)*nbinsrad;
     i2 = i*nbinsrad;
     nbubbles_ax(i) = sum(nbubbles(i1:i2)); 
     if nbubbles_ax(i) > 0 
        bubblevolavg_ax(i) = sum(bubblevolavg(i1:i2).*nbubbles(i1:i2))/nbubbles_ax(i); 
        bubbleAR1avg_ax(i) = sum(bubbleAR1avg(i1:i2).*nbubbles(i1:i2))/nbubbles_ax(i); 
        bubbleCSavg_ax(i) = sum(bubbleCSavg(i1:i2).*nbubbles(i1:i2))/nbubbles_ax(i); 
        bubblecordavg_ax(i) = sum(bubblecordavg(i1:i2).*nbubbles(i1:i2))/nbubbles_ax(i); 
     end
     nbubbles_vel_ax(i) = sum(nbubbles_vel(i1:i2)); 
     if nbubbles_vel_ax(i) > 0 
        bubblevxavg_ax(i) = sum(bubblevxavg(i1:i2).*nbubbles_vel(i1:i2))/nbubbles_vel_ax(i); 
        bubblevyavg_ax(i) = sum(bubblevyavg(i1:i2).*nbubbles_vel(i1:i2))/nbubbles_vel_ax(i); 
     end
end

% aggregate axial bins for radial variation 
nbubbles_rad = zeros(nbinsrad,1); 
bubblevolavg_rad = zeros(nbinsrad,1); 
bubbleAR1avg_rad = zeros(nbinsrad,1); 
bubbleCSavg_rad = zeros(nbinsrad,1); 
bubblecordavg_rad = zeros(nbinsrad,1);
bubblevxavg_rad = zeros(nbinsrad,1);
bubblevyavg_rad = zeros(nbinsrad,1); 
nbubbles_vel_rad = zeros(nbinsrad,1); 

for i=1:nbinsrad 
     i2 = nbinsax*nbinsrad; 
     nbubbles_rad(i) = sum(nbubbles(i:nbinsrad:i2)); 
     if nbubbles_rad(i) > 0 
        bubblevolavg_rad(i) = sum(bubblevolavg(i:nbinsrad:i2).*nbubbles(i:nbinsrad:i2))/nbubbles_rad(i); 
        bubbleAR1avg_rad(i) = sum(bubbleAR1avg(i:nbinsrad:i2).*nbubbles(i:nbinsrad:i2))/nbubbles_rad(i);  
        bubbleCSavg_rad(i) = sum(bubbleCSavg(i:nbinsrad:i2).*nbubbles(i:nbinsrad:i2))/nbubbles_rad(i); 
        bubblecordavg_rad(i) = sum(bubblecordavg(i:nbinsrad:i2).*nbubbles(i:nbinsrad:i2))/nbubbles_rad(i); 
     end
     nbubbles_vel_rad(i) = sum(nbubbles_vel(i:nbinsrad:i2)); 
     if nbubbles_vel_rad(i) > 0 
        bubblevxavg_rad(i) = sum(bubblevxavg(i:nbinsrad:i2).*nbubbles_vel(i:nbinsrad:i2))/nbubbles_vel_rad(i); 
        bubblevyavg_rad(i) = sum(bubblevyavg(i:nbinsrad:i2).*nbubbles_vel(i:nbinsrad:i2))/nbubbles_vel_rad(i); 
     end
end

bubblediaavg = ((4/pi)*bubblevolavg).^(1/2); 
bubblediaavg_ax = ((4/pi)*bubblevolavg_ax).^(1/2); 
bubblediaavg_rad = ((4/pi)*bubblevolavg_rad).^(1/2); 


% 2D histogram 

bubblestats_2D = [bincenter_2Dr, bincenter_2Dy, nbubbles, bubblediaavg, bubbleCSavg, bubblecordavg, bubbleAR1avg, nbubbles_vel, bubblevxavg, bubblevyavg]; 
bubblestats_ax = [bincenter_ax, nbubbles_ax, bubblediaavg_ax, bubbleCSavg_ax, bubblecordavg_ax, bubbleAR1avg_ax, nbubbles_vel_ax, bubblevxavg_ax, bubblevyavg_ax]; 
bubblestats_rad = [bincenter_rad, nbubbles_rad, bubblediaavg_rad, bubbleCSavg_rad, bubblecordavg_rad, bubbleAR1avg_rad, nbubbles_vel_rad, bubblevxavg_rad, bubblevyavg_rad]; 

end



