clear
close all

addpath(genpath('/nobackup/dcarrol2/MATLAB'));

gridDir = '/nobackup/dcarrol2/grid/LLC_270/'
codeDir = '/nobackup/dcarrol2/bin_average/m_files/';
dataDir =  '/nobackup/dcarrol2/bin_average/mat/';
saveDir =  '/nobackup/dcarrol2/bin_average/mat/';

corners = load([dataDir 'cell_corners.mat']);

%%

nx = 270;
ny = nx .* 13;

%%

%input grid
input.XG = readbin([gridDir 'XG.data'],[nx ny]);
input.XG2 = input.XG;
input.XG2(find(input.XG<0)) = input.XG(find(input.XG<0)) + 360;
input.YG = readbin([gridDir 'YG.data'],[nx ny]);
input.dXG = readbin([gridDir 'DXG.data'],[nx ny]);
input.dYG = readbin([gridDir 'DYG.data'],[nx ny]);
input.RAC = readbin([gridDir 'RAC.data'],[nx ny]);
input.hFacC = readbin([gridDir 'hFacC.data'],[nx ny]);
input.AngleCS = readbin([gridDir 'AngleCS.data'],[nx ny]);
input.AngleSN = readbin([gridDir 'AngleSN.data'],[nx ny]);

rSphere = mmax(input.dYG).*4.*nx/2/pi;
rDeg = 2.*pi.*rSphere/360;

%output grid
output.XG = (-180:179)' .* ones(1,180);
output.XG2 = output.XG;
output.XG2(find(output.XG<0)) = output.XG(find(output.XG<0)) + 360;
output.YG = ones(360,1) .* (-90:89);
output.RAC = rDeg^2*cosd(output.YG+.5);

corners.XG2sw(find(corners.XGsw < 0)) = corners.XGsw(find(corners.XGsw < 0)) + 360;
corners.XG2se(find(corners.XGse < 0)) = corners.XGse(find(corners.XGse < 0)) + 360;
corners.XG2ne(find(corners.XGne < 0)) = corners.XGne(find(corners.XGne < 0)) + 360;
corners.XG2nw(find(corners.XGnw < 0)) = corners.XGnw(find(corners.XGnw < 0)) + 360;

%%

% bin averaging template
bin_average = spalloc(length(output.XG(:)),length(input.XG(:)),10*length(input.XG(:)));

inYG = input.YG;
outYG = output.YG;

%%

for i=1:length(output.XG(:))
    
    disp(num2str(output.YG(i)));
    
    %process 90W-90E for input grid on facets 1-2 using XG
    if output.XG(i) >= -90 && output.XG(i) < 90
        
        inXG = input.XG;
        outXG = output.XG;
        
    else %process 90E-270E for input grid on facets 1-2 using XG2
        
        inXG = input.XG2;
        outXG = output.XG2;
        
    end
    
    %find input grid cells that may intersect output grid cell
    ix = find(inXG > (outXG(i)-.34) & inXG <= (outXG(i)+1.34) & ...
        inYG > (outYG(i)-.34) & inYG < (outYG(i)+1.34));
    
    %process the lat/lon, 70S-57N region
    if output.YG(i) >= -70 && output.YG(i) < 57
        
        %convert to rectangular coordinates and define rectangles [X,Y,WIDTH,HEIGHT]
        A = [(inXG(ix)-outXG(i))*rDeg*cosd(outYG(i)+.5) (inYG(ix)-outYG(i))*rDeg ...
            input.dXG(ix)*cosd(outYG(i)+.5)./cosd(inYG(ix)) input.dYG(ix)];
        
        % apply corrections for input grid on facets 4-5
        i45 = find(round(input.AngleCS(ix)) == 0);
        
        if length(i45) > 0
            
            A(i45,2) = A(i45,2) - input.dXG(ix(i45));
            A(i45,3) = input.dYG(ix(i45))*cosd(output.YG(i)+.5)./cosd(input.YG(ix(i45)));
            A(i45,4) = input.dXG(ix(i45));
            
        end
        
        %define output grid rectangle [X,Y,WIDTH,HEIGHT]
        B = [0 0 rDeg*cosd(output.YG(i)+.5) rDeg];
        
        %find intersection areas
        intArea = rectint(A,B);
        
        %compute bin average weight
        if abs(output.RAC(i)-sum(intArea))/output.RAC(i) > 1e-4
            
            error(['output.RAC(' int2str(i) ') differs from sum(intArea)'])
            
        else
            
            intArea=intArea / sum(intArea);
            bin_average(i,ix) = intArea;
            
        end
        
    elseif outYG(i) < - 70 || outYG(i) >= 57 %Arctic polar cap and southern tripolar grid
        
        if(outYG(i) == -90 || outYG(i) == 90)
            
            idx = 10^2;
            idy = 10^2;
        
        else
            
            idx = 0.5 .* 10^3;
            idy = 0.5 .* 10^3;
            
        end
        
        %find input grid cells that may intersect output grid cell
        ix = find(inYG > (outYG(i) - 1.1) & inYG < (outYG(i) + 1.1));
        
        [bx by] = polarstereo_fwd([outYG(i) outYG(i) outYG(i)+1 outYG(i)+1 outYG(i)], ...
            [outXG(i) outXG(i)+1 outXG(i)+1 outXG(i) outXG(i)]);
        
        [x y] = meshgrid([min(bx):idx:max(bx)],[min(by):idy:max(by)]);
        
        in1 = find(inpolygon(x,y,bx,by) == 1); %number of fine grid points in output grid polygon
        
        subX = x(in1);
        subY = y(in1);
        
        in2 = [];
        
        hold on
        
        for i = 1:length(ix)
            
            [cx cy] = polarstereo_fwd([corners.YGsw(ix(i)) corners.YGse(ix(i)) corners.YGne(ix(i)) corners.YGnw(ix(i)) corners.YGsw(ix(i))], ...
                [corners.XGsw(ix(i)) corners.XGse(ix(i)) corners.XGne(ix(i)) corners.XGnw(ix(i)) corners.XGsw(ix(i))]);
            
            in2 = find(inpolygon(subX,subY,cx,cy) == 1);

            intArea(i) = length(in2) / length(in1);
            
        end
        
        bin_average(i,ix) = intArea;
        
        clear intArea
        
    end
    
    clear in1 in2 clear ix
    
end

%%

save([saveDir 'bin_average.mat'],'input','output','bin_average');

%% 
