function [ counterOfIntensity ] = HistogramValues( img )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
counterOfIntensity=zeros(1,256);


for x=1:size(img,1)
    for y=1:size(img,2)
        for intensity=0:255
            if(img(x,y)==intensity)
                counterOfIntensity(intensity+1)=counterOfIntensity(intensity+1)+1;
            end
        end
    end
end


end

