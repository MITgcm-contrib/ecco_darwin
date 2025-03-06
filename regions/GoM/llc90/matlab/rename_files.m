% This code snipet was used to rename files after initially running
% ExtractFields.m with DeltaT=1200 instead 3600.

cd ~/ecco_darwin/GoM
suf1='_20x15.';
suf2='_20x15x47.';
DeltaT=3600;
for fot={'ETAN','apCO2', 'fugfCO2', 'CO2_flux','SIarea','SIheff'}
    dnm=dir([fot{1} '/' fot{1} '*']);
    for t=length(dnm):-1:1
        fin=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fin,suf1);
        ts1=fin((l+length(suf1)):(l+length(suf1))+14);
        ts=(datenum(ts1,'yyyymmddTHHMMSS')-datenum(1992,1,1))*24*60*60/1200;
        dy=ts2dte(ts,DeltaT,1992,1,1,30);
        fout=[dnm(t).folder '/' fot{1} suf1 '.' dy];
        eval(['!mv ' fin ' ' fout])
    end
end
for fot={'THETA', 'DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', ...
         'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', ...
         'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', ...
         'c6', 'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5', 'pH', ...
         'U','V','SALT'}
    dnm=dir([fot{1} '/' fot{1} '*']);
    for t=length(dnm):-1:1
        fin=[dnm(t).folder '/' dnm(t).name];
        l=strfind(fin,suf2);
        ts1=fin((l+length(suf2)):(l+length(suf2))+14);
        ts=(datenum(ts1,'yyyymmddTHHMMSS')-datenum(1992,1,1))*24*60*60/1200;
        dy=ts2dte(ts,DeltaT,1992,1,1,30);
        fout=[dnm(t).folder '/' fot{1} suf1 '.' dy];
        eval(['!mv ' fin ' ' fout])
    end
end
