function [ MaskedArray ] = myRegionGrowAlgorithm( img,T,seed )
%myRegionGrowAlgorithm Implementation of the Region Growing Algorithm
%   Detailed explanation goes here
imgDouble=im2double(img);
[sizeX,sizeY]=size(img);
%Create the Mask and initialize it with zeros
MaskedArray=zeros(sizeX,sizeY);
% Create an empty list where we put the border pixels
BorderPixelsList=[];
%Set the mask to 1 for the seed pixel
MaskedArray(seed(1,1),seed(1,2))=1;
%Add the seed to the borders list
BorderPixelsList=[BorderPixelsList;seed];

%While there is still a border pixel we have to check
while( ~isempty(BorderPixelsList) )
    %Check the validity of the border pixels
    %(borders outside of the image grid)
    validBorderPixels=findValidNeighbors(BorderPixelsList(1,1),BorderPixelsList(1,2),sizeX,sizeY);
    %For each valid neighbor, check if the condition is true to add it to
    %the region
    for currentBorderPixel=1:size(validBorderPixels,1)
        %Another condition could be
        %if(abs(imgDouble(BorderPixelsList(1,1),BorderPixelsList(1,2))-imgDouble(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2)))<=T && img(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2))>=T && MaskedArray(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2))~=1)% less than threshold & not already segmented
        if( img(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2))>=T ...
                && MaskedArray(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2))~=255)  % less than threshold & not already segmented
            MaskedArray(validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2))=255;
            BorderPixelsList=[BorderPixelsList;validBorderPixels(currentBorderPixel,1),validBorderPixels(currentBorderPixel,2)];
        end
    end
    %Delete the first element of the border list, since we just checked it
    BorderPixelsList(1,:)=[];
end


end

