% Look at some output from the global cutout

% Specify location of global cutout, e.g.,
pin='~/Links/Box/Public/GulfGuinea/GlobalCutout/';

% Grid dimensions
nx=69;
ny=66;
nz=50;

% Plot bathymetry
XC=readbin([pin 'grid/XC_69x66'],[69 66]);
YC=readbin([pin 'grid/YC_69x66'],[69 66]);
Depth=readbin([pin 'grid/Depth_69x66'],[69 66]);
Depth(find(Depth==0))=nan;
pcolorcen(XC,YC,-Depth)
axis([min(XC(:)) max(XC(:)) min(YC(:)) max(YC(:))])
colorbar
title('Model bathymetry (m)')
xlabel('Longiude East (^o)')
ylabel('Latitude North (^o)')
