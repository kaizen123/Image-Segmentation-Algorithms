function out=findValidNeighbors(xInput,yInput,sizeX,sizeY)

    out=[];
    if(yInput>1),  out=[out;xInput,yInput-1];   end
    if(xInput>1),  out=[out;xInput-1,yInput];   end
    if(xInput<sizeX),  out=[out;xInput+1,yInput];   end
    if(yInput<sizeY),  out=[out;xInput,yInput+1];   end
    
    
end