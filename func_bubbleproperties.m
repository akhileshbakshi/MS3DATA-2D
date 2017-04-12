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

B(end,:) = []; 
C(end,:) = []; 

% sortrows based on ylocation 
C = sortrows(C,4); 
C(:,2) = []; 

% renumber bubbles 
C(:,1) = linspace(1,length(C(:,1)),length(C(:,1)))';

% output = [bubble#, xmean, ymean, bubbleA, xmin, xmax, ymin, ymax, AR]

end

