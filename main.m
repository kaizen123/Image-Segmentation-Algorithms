clc;
clear all;
close all;

% Read the Girl Face Image(input_image.png)
% and Lena Image(lena512.bmp)
lenaImg=imread('lena512.bmp');
girlImg=imread('input_image.png');

%%
%Task 1.1-Part 1
%Draw the 2 images
fig1=figure;
scrsz = get(0,'ScreenSize');
set(fig1, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
fig1.Name='Task 1.1: Histogram of Girl Face';
%subplot(2,2,1:2);
subplot(2, 2, [1 3]);
imshow(girlImg);
title('Girl Face Image')
IntensitiesRange=0:255;
defaultImgHistValues=HistogramValues(girlImg);
subplot(2,2,2)
plot(IntensitiesRange,defaultImgHistValues,'b');
axis([0 255 0 inf])
xlabel('Intensities');
ylabel('Number of Pixels');
title('2D-Line Plot of Histogram')
subplot(2,2,4)
histo = histogram(girlImg,255,'BinWidth',1);
axis([0 255 0 inf])
xlabel('Intensities');
ylabel('Number of Pixels');
title('Bar Histogram')

%%
%Task 1.1-Part 2
%Plot Histogram Comparison
lenaImgHistValues=HistogramValues(lenaImg);
fig2=figure;
set(fig2, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
fig2.Name='Histogram Comparison';
subplot(3,2,1)
imshow(girlImg);
subplot(3,2,2)
imshow(lenaImg);
subplot(3,2,[3 4])
plot(IntensitiesRange,defaultImgHistValues,'b');
hold on
plot(IntensitiesRange,lenaImgHistValues,'r');
hold off
axis([0 255 0 inf])
legend('Girl Face','Lena')
title('2D-Line Plot of Histogram')
subplot(3,2,[5 6])
%Comparison of Histograms
h1 = histogram(girlImg,255,'BinWidth',1);
hold on
h2 = histogram(lenaImg,255,'BinWidth',1);
hold off
axis([0 255 0 inf])
legend('Girl Face','Lena')
title('Bar Histogram')

%%
%Task 1.2
% Write a short program to threshold the image and try to
%identify a good threshold by trial and error.

fig3=figure;
fig3.Name='Task 1.2: Image with different thresholds';
set(fig3, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
subplot(2,3,1)
imshow(girlImg);
title('Original Image');
subplot(2,3,2)
ThresholdValue=70;
ThresholdedImage=thresholdImage(girlImg,ThresholdValue );
imshow(ThresholdedImage);
title('Thresholded Image-T=70');
ThresholdValue=90;
ThresholdedImage=thresholdImage(girlImg,ThresholdValue );
subplot(2,3,3)
imshow(ThresholdedImage);
title('Thresholded Image-T=90');
ThresholdValue=110;
ThresholdedImage=thresholdImage(girlImg,ThresholdValue );
subplot(2,3,4)
imshow(ThresholdedImage);
title('Thresholded Image-T=110');
ThresholdValue=130;
ThresholdedImage=thresholdImage(girlImg,ThresholdValue );
subplot(2,3,5)
imshow(ThresholdedImage);
title('Thresholded Image-T=130');
ThresholdValue=150;
ThresholdedImage=thresholdImage(girlImg,ThresholdValue );
subplot(2,3,6)
imshow(ThresholdedImage);
title('Thresholded Image-T=150');

%%
%Task B-Part 2: Create a ground truth segmentation.
%You can do this using Matlab (look up the function roipoly),
%or using other painting tools. (15 points)

%This part of the code is used to cut off the ground truth part
%of the given image by the user!
%it is not enabled here,because we get the ground truth image loaded.
%figure,
%BW = roipoly(img);
%figure, imshow(img);
%figure, imshow(BW);
%imwrite(BW,'GroundTruthImage.png');

%%
%Task 1.3
groundTruth=imread('GroundTruthImage.png');
%Change thresholded image to range 0-1
tic;
%Calculate the variables needed to build the ROC Curve
[TPR FPR P N TP FN TN FP]=GetROC(girlImg,groundTruth);
toc

%Calculate the distance between each point of the ROC curve and the
%point (0,1) which is the optimal.
for T=0:255
    distanceFromOptPoint(T+1)=sqrt( ((TPR(T+1)-0)^2)+((FPR(T+1)-1)^2));
end

%Approximate the gradient of each point of the ROC Curve
%by calculating the gradient between each point and the next one on
%the curve
gradientValues=zeros(256,1);
for ind=0:255
    if(ind+1==256 || ind+1==1)
        gradientValues(ind+1)=0;
    else
        gradientValues(ind+1)=TPR((ind+1)+1)-TPR((ind+1))/(FPR((ind+1)+1)-FPR((ind+1)));
    end
end

%Calculate the offset between each point's gradient and
%the correct gradient N/P
optimalGradientVal=N(1)/P(1);
diffFromCorrectGradient=abs( gradientValues-optimalGradientVal );

%Normalize the 2 values
NormDistFromOptimalPoint=(distanceFromOptPoint-min(distanceFromOptPoint))/(max(distanceFromOptPoint)-min(distanceFromOptPoint));
NormDiffFromCorrectGradient=(diffFromCorrectGradient-min(diffFromCorrectGradient))/(max(diffFromCorrectGradient)-min(diffFromCorrectGradient));
%Get the transposed matrix to match with the other array above
NormDiffFromCorrectGradient=NormDiffFromCorrectGradient';

%Calculate a metric that takes into account both the distance
%of each point of the ROC curve from the (0,1) point and also
%the difference between the each point's gradient and the N/P gradient
%The minimum of this metric has the required point
metricValue=0.5*NormDistFromOptimalPoint+0.5*NormDiffFromCorrectGradient;
[MinimumValue,OperatingPoint]=min(metricValue);

%Plot the ROC curve and the operating point
fig4=figure;
fig4.Name='Task 1.3: ROC Curve';
set(fig4, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
plot(FPR,TPR,'b',FPR(OperatingPoint),TPR(OperatingPoint),'r*');
title('ROC Curve')
hold on
plot(0:0.001:1,0:0.001:1,'k');
xlabel('FPR');
ylabel('TPR');

%Plot the thresholded image with T equal to the T of the operating point
fig5=figure;
set(fig5, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
ThresholdedImage=thresholdImage(girlImg,OperatingPoint-1);
subplot(1,2,1);
imshow(girlImg);
title('Original Image');
subplot(1,2,2);
imshow(ThresholdedImage);
str = sprintf('Thresholded Image using the Operating Point: %d',OperatingPoint-1);
title(str);





%%
%Advanced Section
%Region Growing Algorithm
T=95;
% 
% %Get a nx2 vector of selected points from the users
% fig6=figure;
% fig6.Name='Task 1.3: ROC Curve';
% set(fig6, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
% imshow(girlImg);
% [ySel, xSel] = getpts
ySel(1)=257;
xSel(1)=422;

for i=1:size(xSel,1)
    ySel(i)=uint16(ySel(i));
    xSel(i)=uint16(xSel(i));
end
seeds=[xSel,ySel];

NumberOfMaskedImages=size(ySel,1);
%NumberOfWindows=size(ySel,1)+1;

%fig8=figure;
%subplot(3,floor((NumberOfWindows+1)/2),1);
% imshow(girlImg);
% title('Original Image');

%If the user enters more than 1 seed, then we take the union of all of them
%as our final masked image
UnionMaskedImg=zeros(size(girlImg,1),size(girlImg,2));
for i=1:NumberOfMaskedImages
    MaskedImg=zeros(size(girlImg,1),size(girlImg,2));
    MaskedImg=myRegionGrowAlgorithm( girlImg,T,seeds(i,:) );
    MaskedImg=changeRangeOfImage( MaskedImg );
    UnionMaskedImg=UnionMaskedImg | MaskedImg;
end

figure,
subplot(2,1,1);
imshow(girlImg);

hold on;
for i=1:size(xSel)
    plot(ySel(i),xSel(i),'r.','MarkerSize',15);
end
hold off;
title('Original Image');
%imshow(girlImg);

subplot(2,1,2);
imshow(UnionMaskedImg);
title('Final Region Grow Result for T=95');


%Create ROC for Region-Growing
% figure,
% imshow(girlImg);
% [ySel, xSel] = getpts;
% for i=1:size(xSel,1)
%     ySel(i)=uint16(ySel(i));
%     xSel(i)=uint16(xSel(i));
% end
% seeds=[xSel,ySel];
% imshow(girlImg);
% title('Original Image');
fig7=figure;
fig7.Name='Task 2.1: The ROC Comparison';
set(fig7, 'Position', [0 0 scrsz(1,3)/2 scrsz(1,4)]);
subplot(2,1,1);
imshow(girlImg);

NumberOfMaskedImages=size(ySel,1);
NumberOfWindows=size(ySel,1)+1;

RG_TP=zeros(256,1);
RG_FN=zeros(256,1);
RG_TN=zeros(256,1);
RG_FP=zeros(256,1);


%Calculate the region growing algorithm result for thresholds from 0 to 255
for T=0:255
    
    sprintf('Loading Region Growing ROC... %1.1f %%',(T/255)*100)
    UnionMaskedImg=zeros(size(girlImg,1),size(girlImg,2));
    for i=1:NumberOfMaskedImages
        MaskedImg=zeros(size(girlImg,1),size(girlImg,2));
        MaskedImg=myRegionGrowAlgorithm( girlImg,T,seeds(i,:) );
        MaskedImg=changeRangeOfImage( MaskedImg );
        UnionMaskedImg=UnionMaskedImg | MaskedImg;
    end
    
    %Calculate the TP,FN,FP,TN for the Region Growing ROC
    for x=1:size(girlImg,1)
        for y=1:size(girlImg,2)
            if( groundTruth(x,y)==1 && UnionMaskedImg(x,y)==1 )
                RG_TP(T+1)=RG_TP(T+1)+1;
            elseif( groundTruth(x,y)==1 && UnionMaskedImg(x,y)==0 )
                RG_FN(T+1)=RG_FN(T+1)+1;
            elseif(groundTruth(x,y)==0 && UnionMaskedImg(x,y)==0)
                RG_TN(T+1)=RG_TN(T+1)+1;
            elseif(groundTruth(x,y)==0 && UnionMaskedImg(x,y)==1)
                RG_FP(T+1)=RG_FP(T+1)+1;
            end
        end
    end
end

RG_P=RG_TP+RG_FN;
RG_N=RG_FP+RG_TN;

RG_TPR=RG_TP./RG_P;
RG_FPR=RG_FP./RG_N;

%Render the comparison of ROCs
%fig9=figure;
%fig9.Name='Task 2.1: The ROC Comparison';
%set(fig9, 'Position', [0 0 scrsz(1,3) scrsz(1,4)])
subplot(2,1,2);
plot(FPR,TPR,'b');
hold on;
plot(RG_FPR,RG_TPR,'r');
xlabel('FPR & RG_FPR');
ylabel('TPR & RG_TPR');


%Calculate the area for each and print the result
AreaOfRegionGrowingMethod=trapz(sort(RG_FPR),sort(RG_TPR));
AreaOfThresholdingMethod=trapz(sort(FPR),sort(TPR));
sprintf('Area of Region Growing ROC Curve= %1.2f',AreaOfRegionGrowingMethod)
sprintf('Area of Thresholding ROC Curve= %1.2f',AreaOfThresholdingMethod)


%%
%Task 2.2: Mean-shift Algorithm
fig8=figure;
fig8.Name='Task 2.2: Mean-shift Algorithm';
set(fig8, 'Position', [0 0 scrsz(1,3) scrsz(1,4)]);
subplot(3,2,1);
imshow(girlImg);
title('Original Image');
for radius=10:10:50
    %radius=75;
    diameter=2*radius+1;
    %Set the mean Value and the buffer value to a negative big number as 
    %an initial value
    meanVal=-999;
    buffer=-888;
    MeanShiftIntensities=zeros(256,1);
    flagVal=0;

    %For each intensity i from 0-255, calculate the new mean and add it
    % to another array
    for i=0:255
        bi=i;
        %While the meanVal hasn't converged to a specific value
        while (meanVal~=buffer)
            buffer=meanVal;
            if flagVal==1
                bi=buffer;
            end
            nom=0.0;
            denom=0.0;
            %For the window around bi value, calculate the mean value
            for k=bi-floor(diameter/2):1:bi+floor(diameter/2)
                %We check for valid values of intensity k
                if (k>=0 && k<=255)
                    nom=nom+k*defaultImgHistValues(k+1);
                    denom=denom+defaultImgHistValues(k+1);
                end
            end
            %avoid division by 0
            if (nom==0 || denom==0)
                meanVal=0;
            else
                meanVal=round(nom/denom);
            end
            flagVal=1;
            %bi=bi+1;
        end
        flagVal=0;
        MeanShiftIntensities(i+1,1)=uint8(meanVal);
        meanVal=-999;
        buffer=-888;
    end

    %Create the Mean Shift image result,by replacing each pixel's intensity
    %with a new intensity based on the Mean of that intensity in
    %MeanShiftIntensities array
    MeanShiftImg=zeros(size(girlImg,1),size(girlImg,2));
    for x=1:size(girlImg,1)
        for y=1:size(girlImg,2)
            MeanShiftImg(x,y)=uint8(MeanShiftIntensities(girlImg(x,y)+1,1));
        end
    end
    MeanShiftImg=uint8(MeanShiftImg);
    subplot(3,2,1+floor(radius/10));
    imshow(im2double(MeanShiftImg));
    str = sprintf('Radius= %d',radius);
    title(str);
    
end
