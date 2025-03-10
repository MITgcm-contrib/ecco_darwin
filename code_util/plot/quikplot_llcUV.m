function quikplot_llcUV(fldU,fldV);

% Function quikplot_llc(fld)
% plot lat-lon field
%
% INPUTS
% fld  input array of dimension nx*nx*13

if nargin < 2, error('please specify field to plot'); end

[nx ny]=size(fldU);
ny=nx*3;

for i=1:2
if i==1; fld=fldU; else fld=fldV; end
% read face 1, nx*nx*3
f1=fld(:,1:(nx*3));

% read face 2, nx*nx*3
f2=fld(:,(nx*3+1):(nx*6));

% read face 3, nx*nx
f3=fld(:,(nx*6+1):(nx*7));

% read face 4, nx*3*nx
f4=f1';
for f=8:10
    i1=(1:nx)+(f-8)*nx;
    i2=(1:3:(nx*3))+7*nx+f-8;
    f4(i1,:)=fld(:,i2);
end

% read face 5, nx*3*nx
f5=f1';
for f=11:13
    i1=(1:nx)+(f-11)*nx;
    i2=(1:3:(nx*3))+10*nx+f-11;
    f5(i1,:)=fld(:,i2);
end

% plot field
f=nan*ones(4*nx,ny);
f(1:nx,1:ny)=f1;
f((nx+1):(2*nx),1:ny)=f2;
%f(1:nx,(ny+1):(ny+nx/2))=rot90(f3(1:(nx/2),:),1);
%f((2*nx+1):(3*nx),(ny+1):(ny+nx/2))=rot90(f3((nx/2+1):nx,:),3);
f((2*nx+1):(3*nx),1:ny)=rot90(f4,3);
f((3*nx+1):(4*nx),1:ny)=rot90(f5,3);

if i==1; uu=f; else vv=f; end
end
	ue=uu; vn=vv;
	l2=(1:2*nx)+2*nx;
        ue(l2,:)     =  vv(l2,:);
        vn(l2,2:end) = -uu(l2,1:end-1);
        vn(l2,1)     = -uu(l2,1);
	
%quikpcolor(f');
subplot(211)
pcolorcen(ue')
shading flat
set(gca,'xtick',[],'ytick',[])
subplot(212)
pcolorcen(vn')
shading flat
set(gca,'xtick',[],'ytick',[])
