% Look at some output from global cutout
cd ~/Links/Box/Public/GulfGuinea/GlobalCutout
nx=69;
ny=66;
nz=50;

% plot bathymetry
XC=readbin('grid/XC_69x66',[69 66]);
YC=readbin('grid/YC_69x66',[69 66]);
Depth=readbin('grid/Depth_69x66',[69 66]);
Depth(find(Depth==0))=nan;
pcolorcen(XC,YC,-Depth)
axis([mmin(XC) mmax(XC) mmin(YC) mmax(YC)])
colorbar
title('Model bathymetry (m)')
xlabel('Longiude East (^o)')
ylabel('Latitude North (^o)')
