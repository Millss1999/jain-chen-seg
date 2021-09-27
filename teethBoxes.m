%load image
I = imread("teeth_sample.png");
%apply constrast filter
I = adapthisteq(I);
I = rescale(I);
%figure, imshow(I)

%local min
minima = [];
II = [];

%sliding window over image cols - 16 windows of 32px
for i = 1:16 %iterate windows
    windowRowSums = [];
    rowNum = []
    for j =1:472 %iterate each window rows
        rowSum = 0;
        for k = 1:32 %iterate px in each row
            rowSum = rowSum + I(j,((i-1) * k) + k);
        end
        windowRowSums = [windowRowSums; rowSum]; %sum of current row
        rowNum = [rowNum; (length(rowNum) +1)]; % parallel array of row nums
        II = cat(2,windowRowSums,rowNum); %concatenate to get II
    end %end row iteration
    %find local minima for given window
    localMin = islocalmin(II(:,1));
    II = cat(2, II, localMin); %cat to II
    
    %isolate local mins
    isolatedLocalMin = [];
    for row = 1:length(II(:,1))
        if(II(row,3) == 1)
            isolatedLocalMin = [isolatedLocalMin; II(row,:)];
        end
    end
            
    %find optimal minima for window
    y_hat = 210; %manually set y_hat
   
    for iter = 1:length(isolatedLocalMin) %for each local min
        probY = pY(isolatedLocalMin(iter,2), y_hat, 5000);
        maxDepth = max(isolatedLocalMin(:,1));
        probD = pD(isolatedLocalMin(iter,1), 1, maxDepth);
        pYpD = probY*probD;
        isolatedLocalMin(iter,4) = pYpD; %store pYpD
    end
    %find max pYpD
    [maxProb, index] = max(isolatedLocalMin(:,4));
    %append min row to optimal minima list
    minima = [minima; isolatedLocalMin(index, 2)]; 
    
end %end window iteration

windowMidPts = []
for i = 1:16
    windowMidPts = [windowMidPts; ((32 * (i-1)) + 16)];
end
%fit polynomial of (row, mid-window) pts
jawGap =polyfit(windowMidPts, minima, 2);
%figure, imshow(I), hold on, plot(windowMidPts, minima, "r")

close all %end horizontal gap find

%##### START TOOTH 1 #####
%find vertical gaps for two teeth
%local min
minima2 = [];
II2 = [];
figure,imshow(I)
[yHatX, yHatY] = ginput(1); %manually set y_hat

%sliding window over image rows - 8 windows of 59px thickness
for i = 1:8 %iterate windows
    windowColSums = [];
    colNum = [];
    for j =1:70 %iterate each window cols -> NOTE: limiter = 70 for purpose of assignment(only 2 teeth)
        colSum = 0;
        for k = 1:59 %iterate px in each col
            colSum = colSum + I((((i-1) * k) + k), j);
        end
        windowColSums = [windowColSums; colSum]; %sum of current col
        colNum = [colNum; (length(colNum) +1)]; % parallel array of col nums
        II2 = cat(2,windowColSums,colNum); %concatenate to get II
    end %end col iteration
    %find local minima for given window
    localMin = islocalmin(II2(:,1));
    II2 = cat(2, II2, localMin); %cat to II
    
    %isolate local mins
    isolatedLocalMin = [];
    for col = 1:length(II2(:,1))
        if(II2(col,3) == 1)
            isolatedLocalMin = [isolatedLocalMin; II2(col,:)];
        end
    end
            
    %find optimal minima for window  
   
    for iter = 1:length(isolatedLocalMin) %for each local min
        probY = pY(isolatedLocalMin(iter,2), yHatX, 5000);
        maxDepth = max(isolatedLocalMin(:,1));
        probD = pD(isolatedLocalMin(iter,1), 1, maxDepth);
        pYpD = probY*probD;
        isolatedLocalMin(iter,4) = pYpD; %store pYpD
    end
    %find max pYpD
    [maxProb, index] = max(isolatedLocalMin(:,4));
    %append min row to optimal minima list
    minima2 = [minima2; isolatedLocalMin(index, 2)]; 
    
end %end window iteration

windowMidPts2 = [];
for i = 1:8
    windowMidPts2 = [windowMidPts2; ((59 * (i-1)) + 8)];
end
%fit polynomial of (row, mid-window) pts
p2 =polyfit(minima2, windowMidPts2, 2);
%figure, imshow(I), hold on, plot(minima2, windowMidPts, "r")
%find end side for tooth 1
%local min
minima2_2 = [];
II2_2 = [];
figure,imshow(I)
[yHatX, yHatY] = ginput(1); %manually set y_hat

%sliding window over image rows - 8 windows of 59px thickness
for i = 1:8 %iterate windows
    windowColSums = [];
    colNum = [];
    for j =70:140 %iterate each window cols -> NOTE: limiter = 70 for purpose of assignment(only 2 teeth)
        colSum = 0;
        for k = 1:59 %iterate px in each col
            colSum = colSum + I((((i-1) * k) + k), j);
        end
        windowColSums = [windowColSums; colSum]; %sum of current col
        colNum = [colNum; (length(colNum) +1 + 70)]; % parallel array of col nums
        II2_2 = cat(2,windowColSums,colNum); %concatenate to get II
    end %end col iteration
    %find local minima for given window
    localMin = islocalmin(II2_2(:,1));
    II2_2 = cat(2, II2_2, localMin); %cat to II
    
    %isolate local mins
    isolatedLocalMin = [];
    for col = 1:length(II2_2(:,1))
        if(II2_2(col,3) == 1)
            isolatedLocalMin = [isolatedLocalMin; II2_2(col,:)];
        end
    end
            
    %find optimal minima for window  
   
    for iter = 1:length(isolatedLocalMin) %for each local min
        probY = pY(isolatedLocalMin(iter,2), yHatX, 5000);
        maxDepth = max(isolatedLocalMin(:,1));
        probD = pD(isolatedLocalMin(iter,1), 1, maxDepth);
        pYpD = probY*probD;
        isolatedLocalMin(iter,4) = pYpD; %store pYpD
    end
    %find max pYpD
    [maxProb, index] = max(isolatedLocalMin(:,4));
    %append min row to optimal minima list
    minima2_2 = [minima2_2; isolatedLocalMin(index, 2)]; 
    
end %end window iteration

windowMidPts2_2 = [];
for i = 1:8
    windowMidPts2_2 = [windowMidPts2_2; (70+ ((59 * (i-1)) + 8))];
end
%fit polynomial of (row, mid-window) pts
p2_2 =polyfit(minima2_2, windowMidPts2_2, 2);
% figure, imshow(I), hold on, plot(minima2_2, windowMidPts2_2, "y"), plot(minima2, windowMidPts2, "y"), plot(windowMidPts, minima, "r")
figure, imshow(I), hold on, plot(minima2_2, windowMidPts2_2, "y"), plot(minima2, windowMidPts2, "y"), plot(1:472, polyval(jawGap, 1:472), "r")

close all %end vertical gap find T1

%##### START T2 Gap find (using t1 gap as first border)
%find end side for tooth 2
%local min
minima3_2 = [];
II3_2 = [];
figure,imshow(I)
[yHatX, yHatY] = ginput(1); %manually set y_hat

%sliding window over image rows - 8 windows of 59px thickness
for i = 1:8 %iterate windows
    windowColSums = [];
    colNum = [];
    for j =140:210 %iterate each window cols -> NOTE: limiter = 70 for purpose of assignment(only 2 teeth)
        colSum = 0;
        for k = 1:59 %iterate px in each col
            colSum = colSum + I((((i-1) * k) + k), j);
        end
        windowColSums = [windowColSums; colSum]; %sum of current col
        colNum = [colNum; (length(colNum) +1 + 140)]; % parallel array of col nums
        II3_2 = cat(2,windowColSums,colNum); %concatenate to get II
    end %end col iteration
    %find local minima for given window
    localMin = islocalmin(II3_2(:,1));
    II3_2 = cat(2, II3_2, localMin); %cat to II
    
    %isolate local mins
    isolatedLocalMin = [];
    for col = 1:length(II3_2(:,1))
        if(II3_2(col,3) == 1)
            isolatedLocalMin = [isolatedLocalMin; II3_2(col,:)];
        end
    end
            
    %find optimal minima for window  
   
    for iter = 1:length(isolatedLocalMin) %for each local min
        probY = pY(isolatedLocalMin(iter,2), yHatX, 5000);
        maxDepth = max(isolatedLocalMin(:,1));
        probD = pD(isolatedLocalMin(iter,1), 1, maxDepth);
        pYpD = probY*probD;
        isolatedLocalMin(iter,4) = pYpD; %store pYpD
    end
    %find max pYpD
    [maxProb, index] = max(isolatedLocalMin(:,4));
    %append min row to optimal minima list
    minima3_2 = [minima3_2; isolatedLocalMin(index, 2)]; 
    
end %end window iteration

windowMidPts3_2 = [];
for i = 1:8
    windowMidPts3_2 = [windowMidPts3_2; (70+ ((59 * (i-1)) + 8))];
end
%fit polynomial of (row, mid-window) pts
p3_2 =polyfit(minima3_2, windowMidPts3_2, 2);
% figure, imshow(I), hold on, plot(minima2_2, windowMidPts2_2, "y"),
figure, imshow(I), hold on, plot(minima2, windowMidPts2, "y"), plot(windowMidPts, minima, "r"), plot(minima2_2, windowMidPts2_2, "y"), plot(minima3_2, windowMidPts3_2, "y")
%figure, imshow(I), hold on, plot(minima2_2, windowMidPts2_2, "y"), plot(minima2, windowMidPts2, "y"), plot(1:472, polyval(jawGap, 1:472), "r"), plot(minima3_2, windowMidPts3_2, "y")

%close all %end vertical gap find T1

