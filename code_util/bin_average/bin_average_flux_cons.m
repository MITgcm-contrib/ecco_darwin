clear
close all;

gridDir = '/Users/carrolld/Documents/research/carbon/simulations/grid/LLC_270/';

codeDir = '/Users/carrolld/Documents/research/carbon/m_files/bin_average/flux_conserving/from_dustin/offline/1x1_deg/';
saveDir = '/Users/carrolld/Documents/research/carbon/mat/bin_average/flux_conserving/';

corners = load([codeDir 'cell_corners.mat']);

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
output.frac = output.XG .* 0;

corners.XG2ne = corners.XGne;
corners.XG2nw = corners.XGnw;
corners.XG2se = corners.XGse;
corners.XG2sw = corners.XGsw;

corners.XG2ne(find(corners.XGne<0)) = corners.XGne(find(corners.XGne<0)) + 360;
corners.XG2nw(find(corners.XGnw<0)) = corners.XGnw(find(corners.XGnw<0)) + 360;
corners.XG2se(find(corners.XGse<0)) = corners.XGse(find(corners.XGse<0)) + 360;
corners.XG2sw(find(corners.XGsw<0)) = corners.XGsw(find(corners.XGsw<0)) + 360;

radius = 6378137.0;
eccentricity = 0.08181919;

%%

% bin averaging template
bin_average_cons = spalloc(length(output.XG(:)),length(input.XG(:)),10*length(input.XG(:)));

inYG = input.YG;
outYG = output.YG;

%%

tic

for i = 1:1:length(output.YG(:))
 
    if(output.YG(i) >= 0)
        
        latTrue = 70;
        
    else
        
        latTrue = -70;
        
    end
    
    %process 90W-90E for input grid on facets 1-2 using XG
    if output.XG(i) >= -90 && output.XG(i) < 90
        
        inXG = input.XG;
        outXG = output.XG;
        
        cXGne = corners.XGne;
        cXGnw = corners.XGnw;
        cXGse = corners.XGse;
        cXGsw = corners.XGsw;
        
    else %process 90E-270E for input grid on facets 1-2 using XG2
        
        inXG = input.XG2;
        outXG = output.XG2;
        
        cXGne = corners.XG2ne;
        cXGnw = corners.XG2nw;
        cXGse = corners.XG2se;
        cXGsw = corners.XG2sw;
        
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
        
        intArea = intArea / sum(intArea);
        bin_average_cons(i,ix) = intArea;
        
    elseif outYG(i) < - 70 || outYG(i) >= 57 %Arctic polar cap and southern tripolar grid
        
        clear ix
        
        %if near pole, use finer spacing for triangle grid geometry
        if(outYG(i) <= -88 || outYG(i) >= 88)
            
            idx = 25;
            idy = 25;
            
        else
            
            idx = 0.5 .* 10^3;
            idy = 0.5 .* 10^3;
            
        end
        
        %find input grid cells within latitudinal band
        ix = find(inYG >= (outYG(i) - 1.5) & inYG <= (outYG(i) + 1.5));
        
        %convert to x,y coordinates
        [bx by] = polarstereo_fwd([outYG(i) outYG(i) outYG(i)+1 outYG(i)+1 outYG(i)], ...
            [outXG(i) outXG(i)+1 outXG(i)+1 outXG(i) outXG(i)], ...
            radius,eccentricity,latTrue);
        
        %hold on
        %scatter(bx,by,100,'r','filled');
        %line(bx,by);
        
        %create fine grid
        [x y] = meshgrid([min(bx)-20.*idx:idx:max(bx)+20.*idx],[min(by)-20.*idy:idy:max(by)+20.*idy]);
        
        in1 = find(inpolygon(x,y,bx,by) == 1); %find number of fine grid points in output grid polygon
        
        subX = x(in1);
        subY = y(in1);
        
        dupIndex = [];
        
        for j = 1:length(ix)
            
            %convert each input grid cell to x,y
            [cx cy] = polarstereo_fwd([corners.YGsw(ix(j)) corners.YGse(ix(j)) corners.YGne(ix(j)) corners.YGnw(ix(j)) corners.YGsw(ix(j))], ...
                [cXGsw(ix(j)) cXGse(ix(j)) cXGne(ix(j)) cXGnw(ix(j)) cXGsw(ix(j))], ...
                radius,eccentricity,latTrue);
            
            %line(cx,cy,'Color','k');
            %drawnow
            
            %find number of fine grid cells bounded by input grid cell
            in2 = find(inpolygon(subX,subY,cx,cy) == 1);
            
            %remove duplicates
            in3 = ismember(in2,dupIndex);
            in2(in3) = [];
            dupIndex = [dupIndex; in2];
            
            %fractional area of output grid cell
            intArea(j) = length(in2) / length(subX(:));
            
            %scatter(subX(in2),subY(in2),100,'b','filled');
            
            clear in2 in3
            
        end
        
        bin_average_cons(i,ix) = intArea;
        
    end
    
    output.frac(i) = nansum(intArea);
    
    disp(['lon: ' num2str(output.XG(i)) ', lat: ' num2str(output.YG(i)) ', latTrue: ' ...
        num2str(latTrue) ', intArea: ' num2str(sum(intArea))]);
    
    if(isnan(intArea) | (nansum(intArea) == 0))
        
        %error('no input grid cells found!');
        
    end
    
    if abs(1-nansum(intArea))/1 > 1e-4
        
        disp(['output.RAC(' int2str(i) ') differs from sum(intArea)'])
        
    end
    
    clear intArea in1
    
end

toc

%%

save([saveDir 'bin_average_LLC_270_to_1x1_deg_test.mat'],'input','output','bin_average_cons','-v7.3');

%%
%test bin-averaging code

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

clf
quikplot_llc(input.hFacC);
colorbar('horiz');
caxis([0 1]);

mask = reshape(bin_average_cons*input.hFacC(:),[360 180]);

hFig2 = figure(2);
set(hFig2,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

mypcolor(mask')
colorbar('horiz');
caxis([0 1]);

%%
