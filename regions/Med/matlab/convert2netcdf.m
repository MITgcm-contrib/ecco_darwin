% cd ~/Links/Box/Public/ECCO_Darwin/Med/grid
nx=144;
ny=72;
nz=47;
XC=readbin('XC_144x72',[nx ny]);
YC=readbin('YC_144x72',[nx ny]);
RC=readbin('RC_47',nz);
hFacC=readbin('hFacC_144x72x47',[nx ny nz]);

nccreate('hFacC.nc','XC','Dimensions',{'x',nx,'y',ny},'Datatype','single','Format','classic')
ncwrite('hFacC.nc','XC',XC)
ncwriteatt('hFacC.nc','XC','description', ...
           'Longitude of grid cell center (degrees East)');

nccreate('hFacC.nc','YC','Dimensions',{'x',nx,'y',ny},'Datatype','single','Format','classic')
ncwrite('hFacC.nc','YC',YC)
ncwriteatt('hFacC.nc','YC','description', ...
           'Latitude of grid cell center (degrees North)');

nccreate('hFacC.nc','RC','Dimensions',{'z',nz},'Datatype','single','Format','classic')
ncwrite('hFacC.nc','RC',RC)
ncwriteatt('hFacC.nc','RC','description', ...
           'Depth of grid cell center (meters)');

nccreate('hFacC.nc','hFacC','Dimensions',{'x',nx,'y',ny,'z',nz},'Datatype','single','Format','classic')
ncwrite('hFacC.nc','hFacC',hFacC)
ncwriteatt('hFacC.nc','hFacC','description', ...
           'Landmask (hFacC=0 is land; 0<hFacC<1 is partially wet grid cell; hFacC=1 is fully wet grid cell')
