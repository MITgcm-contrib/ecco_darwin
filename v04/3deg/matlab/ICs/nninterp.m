function [newGridIndex oldGridIndex] = nninterp(oldGrid,newGrid,zLevel,dl)

%reshape into vectors
oldXC = oldGrid.XC(:);
newXC = newGrid.XC(:);

%if near dateline, then use unwrapped grid
oldXC2 = oldXC;
ix1 = find(oldXC2 < 0);
oldXC2(ix1) = oldXC2(ix1) + 360;

newXC2 = newXC;
ix2 = find(newXC2 < 0);
newXC2(ix2) = newXC2(ix2) + 360;

oldYC = oldGrid.YC(:);
newYC = newGrid.YC(:);

oldHFacC = oldGrid.hFacC(:,:,zLevel);
oldHFacC = oldHFacC(:);

newHFacC = newGrid.hFacC(:,:,zLevel);
newHFacC = newHFacC(:);

interpField = newHFacC*0; %placeholder for interpolated field

c = 1;

for i = 1:length(newXC) %loop through new grid
    
    if(newHFacC(i) > 0) %if wet cell then interpolate
        
        newLon = newXC(i);
        newLat = newYC(i);
        
        xi = find(oldXC >= newXC(i)-dl/2 & oldXC <= newXC(i)+dl/2 & oldYC >= newYC(i)-dl/2 & oldYC <= newYC(i)+dl/2);
        
        nnDistance = dl*(111*1000); %initial distance set to be largest
        
        nxi = [];
        
        %minimize distances within bounding box
        for j = 1:length(xi)
            
            oldLon = oldXC(xi(j));
            oldLat = oldYC(xi(j));
            
            dLat = newLat-oldLat;
            dLon = (newLon-oldLon)*cosd(oldLat);
            
            distance = sqrt(dLon.^2 + dLat.^2);
            
            %if smallest distance, save indices
            if (distance < nnDistance && oldHFacC(xi(j)) ~= 0)
                
                nnDistance = distance;
                nxi = j;
                
            end
            
        end
        
        if(~isempty(nxi)) %found wet nn within bounding box
            
            newGridIndex(c) = i; %index in new grid
            oldGridIndex(c) = xi(nxi); %index in old grid
            
        else %did not find wet nn, post process these
            
            newGridIndex(c) = i; %index in new grid
            oldGridIndex(c) = nan; %index in old grid
            
        end
        
        c = c+1;
        
    end
    
end

end
