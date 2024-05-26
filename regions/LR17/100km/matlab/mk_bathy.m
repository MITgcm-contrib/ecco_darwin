[x y z]=grdread2('MNT_PC20m_HOMONIM_WGS84_NM_ZNEG_V2.0.grd');
z(find(isnan(z)))=mmax(z);
figure(1)
pcolorcen(x,y,z);
colorbar

[NY NX]=size(z);
nx=floor((NX-4)/4);
ny=floor(NY/4);

lon=zeros(nx,1);
lat=zeros(ny,1);
bathy=zeros(nx,ny);

for i=1:nx
    lon(i)=mean(x(((i-1)*4+1):(i*4)));
end
for j=1:ny
    lat(j)=mean(y(((j-1)*4+1):(j*4)));
end
for i=1:nx
    for j=1:ny
        ix=((i-1)*4+1):(i*4);
        jx=((j-1)*4+1):(j*4);
        bathy(i,j)=mmean(z(jx,ix));
    end
end

figure(2)
pcolorcen(lon,lat,bathy');
colorbar

writebin('bathy_LR17_936x875',bathy)
