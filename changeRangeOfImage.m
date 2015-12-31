function [ out_img ] = changeRangeOfImage( img )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
out_img=zeros(size(img,1),size(img,2));
for x=1:size(img,1)
    for y=1:size(img,2)
        %If white=255 turn it to 1
        if(img(x,y)==255)
            out_img(x,y)=1;
        else
            out_img(x,y)=0;
        end
    end
end

end

