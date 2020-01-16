function [lon lat bin_average] = compute_bin_average(XC,YC,RAC,LonIncr,LatIncr,saveDir,saveFilename)

XC = XC(:);
YC = YC(:);
RAC = RAC(:);

%% 
%Abhishek code for grid area

% specify latitude and longitude increment
%LonIncr = 4;
%LatIncr = 5;

% calculate lat-long boundaries
lon =-180+LonIncr/2:LonIncr:180;
lat =-90+LatIncr/2:LatIncr:90;

LatIncr = nanmean(diff(lat));
LonIncr = nanmean(diff(lon));

[xx yy] = meshgrid(lon,lat);

nlat = 180/LatIncr;
nlon = 360/LonIncr;

% Calculate gridwise areas
r=6375.*1000;        
fjep=0.5*(nlat+1);
dlat=pi/(nlat-1);
dd=2*r^2*2*pi/nlon*sin(0.5*dlat);

for j=2:nlat-1
    
    dxyp(j) = dd*cos(dlat*(j-fjep));

end

dxyp(1)=2*r^2*2*pi*sin(0.25*dlat)*cos(0.25*(2*pi-dlat))/nlon;
dxyp(nlat)=dxyp(1);

AREA = repmat(dxyp,nlon,1);  % final output

%% 

% define edges of output grid
lat1 = lat - LatIncr/2;
lat1(1) = -90;

lat2 = lat + LatIncr/2;
lat2(end) = 90;

lon1 = lon - LonIncr/2;
lon2 = lon + LonIncr/2;

[LAT1 LON1] = meshgrid(lat1,lon1);
[LAT2 LON2] = meshgrid(lat2,lon2);

%% 

% put XC in same range as LON1 and LON2
ix = find(XC < (min(LON1(:))));

if length(ix) > 0
    
    XC(ix) = XC(ix) + 360;
    
end

clear ix

ix = find(XC >= (max(LON2(:))));

if length(ix) > 0
    
    XC(ix) = XC(ix) - 360;
    
end

%% 

% Compute bin-averaging template
LON1v = LON1(:);
LAT1v = LAT1(:);

LON2v = LON2(:);
LAT2v = LAT2(:);

XCv = XC(:);
YCv = YC(:);

RACv = RAC(:);
AREAv = AREA(:);

bin_average = spalloc(length(LON1v),length(XCv),length(XCv));

for i=1:length(LON1v)
    
    ix = find(XCv >= LON1v(i) & XCv < LON2v(i) & YCv >= LAT1v(i) & YCv < LAT2v(i));
    
    if length(ix) > 0
        
        %bin_average(i,ix) = RACv(ix) / AREAv(i); %flux-conserving bin average
        bin_average(i,ix) = 1/length(ix); %normal bin average
        
    end
    
    disp(num2str(i));
    
end

save([saveDir saveFilename],'lon','lat','AREA','bin_average');

end