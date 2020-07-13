clear

nx=46;ny=68;
mask=zeros([nx ny]);
relaxMaskAFile='Mac_rbcs_mask.bin';

%from data.obcs
%  OB_Jsouth =   46*1,
%  OB_Iwest  =   68*1,
  seaiceSpongeThickness = 5,
  Arelaxobcsinner = 432000.0,
  Arelaxobcsbound = 86400.0,
%from data.rbcs  
   tauRelaxA=432000.,

dt = (Arelaxobcsinner-Arelaxobcsbound)/(seaiceSpongeThickness-1);
Teff = Arelaxobcsbound + ([1:seaiceSpongeThickness]-1)*dt;
mm = tauRelaxA./Teff;

%West sponge
for i=1:seaiceSpongeThickness
	mask(i,:)=mm(i);
end
%South sponge
for i=1:seaiceSpongeThickness
	mask(i:end,i)=mm(i);
end

writebin(relaxMaskAFile,mask)
