clear, close all

% Mac_Delta270 dimensions and grid file location
nx = 46;
ny = 68;
gDir = '/Users/carrolld/Documents/research/mackenzie/grid/LLC_270/';

% llc270 dimensions and tile file location
NX=270;
tDir = gDir;

% input grid and landmask
XGsw = readbin([gDir 'XG.data'],[nx ny]);
YGsw = readbin([gDir 'YG.data'],[nx ny]);

% read tile00?.mitgdrid files
tile{1}.XG=readbin([tDir 'tile001.mitgrid'],[(NX+1) (NX*3+1)],1,'real*8',5);
tile{1}.YG=readbin([tDir 'tile001.mitgrid'],[(NX+1) (NX*3+1)],1,'real*8',6);
tile{2}.XG=readbin([tDir 'tile002.mitgrid'],[(NX+1) (NX*3+1)],1,'real*8',5);
tile{2}.YG=readbin([tDir 'tile002.mitgrid'],[(NX+1) (NX*3+1)],1,'real*8',6);
tile{3}.XG=readbin([tDir 'tile003.mitgrid'],[(NX+1) (NX+1)],1,'real*8',5);
tile{3}.YG=readbin([tDir 'tile003.mitgrid'],[(NX+1) (NX+1)],1,'real*8',6);
tile{4}.XG=readbin([tDir 'tile004.mitgrid'],[(NX*3+1) (NX+1)],1,'real*8',5);
tile{4}.YG=readbin([tDir 'tile004.mitgrid'],[(NX*3+1) (NX+1)],1,'real*8',6);
tile{5}.XG=readbin([tDir 'tile005.mitgrid'],[(NX*3+1) (NX+1)],1,'real*8',5);
tile{5}.YG=readbin([tDir 'tile005.mitgrid'],[(NX*3+1) (NX+1)],1,'real*8',6);

% find the remaining 3 corners
XGse=XGsw; XGnw=XGsw; XGne=XGsw; YGse=YGsw; YGnw=YGsw; YGne=XGsw;
for i=1:length(XGsw(:)), disp(i)
    for t=1:5
        [I J] = find( abs(XGsw(i)-tile{t}.XG(1:end-1,1:end-1))<1e-4 & ...
                      abs(YGsw(i)-tile{t}.YG(1:end-1,1:end-1))<1e-4 );
        if length(I)>0
            break
        end
    end
    
    %    if length(I)~=1
    %        disp([i t])
    %        disp([I J])
    %        break
    %    end
    
    I=I(1); J=J(1);
    XGse(i)=tile{t}.XG(I+1,J);
    XGnw(i)=tile{t}.XG(I,J+1);
    XGne(i)=tile{t}.XG(I+1,J+1);
    YGse(i)=tile{t}.YG(I+1,J);
    YGnw(i)=tile{t}.YG(I,J+1);
    YGne(i)=tile{t}.YG(I+1,J+1);
    
end

save cell_corners X* Y*

% look @ fields
figure(1), clf, pcolorcen(XGse'-XGsw'), colorbar('horiz'), title('XGse-XGsw')
figure(2), clf, pcolorcen(XGne'-XGnw'), colorbar('horiz'), title('XGne-XGnw')
figure(3), clf, pcolorcen(YGnw'-YGsw'), colorbar('horiz'), title('YGnw-YGsw')
figure(4), clf, pcolorcen(YGne'-YGse'), colorbar('horiz'), title('YGne-YGse')
