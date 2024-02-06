function [] = model_obs_compare_profiles(loadData,modelDir)

fileID = fopen('model_obs_compare_profiles_log.txt','w');

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';

dataDir = '/nobackup/dcarrol2/evaluation/mat/observations/GLODAP/';
codeDir = '/nobackup/dcarrol2/evaluation/m_files/model_obs_compare/GLODAP/';
saveDir = '/nobackup/dcarrol2/evaluation/mat/model_obs_compare/GLODAP/';

%%

numFaces = 13;
nx = 270;
ny = 270 .* numFaces;
nz = 50;

XC = readbin([gridDir 'XC.data'],[nx ny],1,'real*4');
YC = readbin([gridDir 'YC.data'],[nx ny],1,'real*4');
hFacC = readbin([gridDir 'hFacC.data'],[nx ny nz],1,'real*4');

landMask = squeeze(hFacC(:,:,1));

%%

if loadData
    
    load([dataDir 'GLODAPv2_2023_profiles.mat']);
    
    time = profiles.time;
    
    fileList = dir([modelDir 'DIC.0*.*data']);
    
    for i = 1:length(fileList)
        
        fileName = fileList(i).name;
        
        temp = datenum(ts2dte(str2num(fileName(end-14:end-5)),1200,1992));

        modelTimeLowerBound(i) = addtodate(double(temp), -1, 'month');
        
        temp = addtodate(double(temp), -1, 'day');
        
        [year month day hour minute seconds] = datevec(temp);
        
        modelTimeUpperBound(i) = datenum(year,month,eomday(year,month),23,59,59);
        
    end
    
    disp(['dumps: ' num2str(length(modelTimeUpperBound))]);
    
    %%
    
    c = 1;
    
    for i = 1:length(modelTimeLowerBound) %for all model output
        
        ti = find(time >= modelTimeLowerBound(i) & time <= modelTimeUpperBound(i)); %find all profiles within model time period
        
        if (~isempty(ti))
            
            lon = profiles.lon(ti);
            lat = profiles.lat(ti);
            
            for j = 1:length(lon) %find corresponding model indices for obs
                
                dist = (((XC - lon(j)).^2) * cosd(lat(j))) + ((YC - lat(j)).^2);
                
                md = find(dist == nanmin(nanmin(dist))); %find index of closest grid cell
                
                if(length(md) > 1) %if profile is equal distance to two grid cells, use first one
                    
                    im(j) = md(1);
                    
                else
                    
                    im(j) = md;
                    
                end
                
            end
            
        else %no profiles exist for model time period
            
            im = [];
            
        end
        
        if (~isempty(im)) %if model index was found then store it
            
            xi = find(landMask(im) == 1);
            
            rawObs{c}.time = time(ti(xi));
            rawObs{c}.lon = profiles.lon(ti(xi));
            rawObs{c}.lat = profiles.lat(ti(xi));
            rawObs{c}.depth = profiles.depth(ti(xi),:);
            rawObs{c}.type = profiles.type(ti(xi));
            
            rawObs{c}.numObs = length(im(xi));
            rawObs{c}.gridIndex = im(xi);
            rawObs{c}.fileIndex = i;
            
            rawObs{c}.THETA = profiles.THETA(ti(xi),:);
            rawObs{c}.SAL = profiles.SAL(ti(xi),:);
            rawObs{c}.NO3 = profiles.NO3(ti(xi),:);
            rawObs{c}.NO2 = profiles.NO2(ti(xi),:);
            rawObs{c}.PO4 = profiles.PO4(ti(xi),:);
            rawObs{c}.SiO2 = profiles.SiO4(ti(xi),:);
            rawObs{c}.O2 = profiles.O2(ti(xi),:);
            
            rawObs{c}.DOC = profiles.DOC(ti(xi),:);
            rawObs{c}.DON = profiles.DON(ti(xi),:);
            
            rawObs{c}.DIC = profiles.TCO2(ti(xi),:);
            rawObs{c}.DIC_pH_25C = profiles. DIC_pH_25C(ti(xi),:);
            rawObs{c}.DIC_pH_insitu = profiles.DIC_pH_insitu(ti(xi),:);
            
            rawObs{c}.ALK = profiles.TALK(ti(xi),:);
            
            rawObs{c}.PH = profiles.pH_insitu(ti(xi),:);
            
            rawObs{c}.CHLA = profiles.CHLA(ti(xi),:);
            
            c = c + 1;
            
        end
        
        clear im
        
        disp(num2str(i));
        
    end
    
    save([saveDir 'profiles_model_grid.mat'],'rawObs','modelTimeUpperBound','modelTimeLowerBound','-v7.3');
    
else
    
    load([saveDir 'profiles_model_grid.mat']);
    
end

%%
%model output

files.THETA = dir([modelDir 'THETA.0*.*data']);
files.SAL = dir([modelDir 'SALTanom.0*.*data']);

files.NO3 = dir([modelDir 'NO3.0*.*data']);
files.NO2 = dir([modelDir 'NO2.0*.*data']);
files.PO4 = dir([modelDir 'PO4.0*.*data']);
files.SiO2 = dir([modelDir 'SiO2.0*.*data']);
files.O2 = dir([modelDir 'O2.0*.*data']);
files.DOC = dir([modelDir 'DOC.0*.*data']);
files.DON = dir([modelDir 'DON.0*.*data']);
files.DIC = dir([modelDir 'DIC.0*.*data']);
files.ALK = dir([modelDir 'ALK.0*.*data']);
files.PH = dir([modelDir 'pH.0*.*data']);
files.CHL1 = dir([modelDir 'Chl1.0*.*data']);
files.CHL2 = dir([modelDir 'Chl2.0*.*data']);
files.CHL3 = dir([modelDir 'Chl3.0*.*data']);
files.CHL4 = dir([modelDir 'Chl4.0*.*data']);
files.CHL5 = dir([modelDir 'Chl5.0*.*data']);

%%

cTHETA = 1;
cSAL = 1;
cNO3 = 1;
cNO2 = 1;
cPO4 = 1;
cSiO2 = 1;
cO2 = 1;
cDOC = 1;
cDON = 1;
cDIC = 1;
cALK = 1;
cPH = 1;
cCHLA = 1;

for i = 1:length(rawObs)
    
    for j = 1:rawObs{i}.numObs %get all observational profiles
        
        clear ziTHETA ziSAL ziNO3 ziNO2 ziPO4 ziSiO2 siO2 ziDOC ziDON ziDIC ziDIC_pH_25C zi_pH_insitu ziALK ziPH ziCHLA
        %find where profiles exist
        ziTHETA = find(~isnan(rawObs{i}.THETA(j,:)));
        ziSAL= find(~isnan(rawObs{i}.SAL(j,:))); 
        ziNO3 = find(~isnan(rawObs{i}.NO3(j,:)));
        ziNO2 = find(~isnan(rawObs{i}.NO2(j,:)));
        ziPO4 = find(~isnan(rawObs{i}.PO4(j,:)));
        ziSiO2 = find(~isnan(rawObs{i}.SiO2(j,:)));
        ziO2 = find(~isnan(rawObs{i}.O2(j,:)));
        ziDOC = find(~isnan(rawObs{i}.DOC(j,:)));
        ziDON = find(~isnan(rawObs{i}.DON(j,:)));
        ziDIC = find(~isnan(rawObs{i}.DIC(j,:)));
        
        %ziDIC_pH_25C = find(~isnan(rawObs{i}.DIC_pH_25C(j,:)));
        %ziDIC_pH_insitu = find(~isnan(rawObs{i}.DIC_pH_insitu(j,:)));
        
        ziALK = find(~isnan(rawObs{i}.ALK(j,:)));
        ziPH = find(~isnan(rawObs{i}.PH(j,:)));
        ziCHLA = find(~isnan(rawObs{i}.CHLA(j,:)));
        
        %%
        %THETA
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziTHETA,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziTHETA(hi) = [];
        
        if (~isempty(ziTHETA))
            
            THETA_obs = rawObs{i}.THETA(j,ziTHETA)';
            
            THETA_path = [files.THETA(rawObs{i}.fileIndex).folder '/' files.THETA(rawObs{i}.fileIndex).name];
            THETA_model = read_ijk(THETA_path,rawObs{i}.gridIndex(j),1,ziTHETA,nx*ny,1,nz);
            
            observations.THETA{cTHETA}.data = THETA_obs;
            observations.THETA{cTHETA}.time = (ones(1,length(THETA_obs)) .* rawObs{i}.time(j))';
            observations.THETA{cTHETA}.depth = rawObs{i}.depth(j,ziTHETA)';

	    %disp(num2str(observations.THETA{cTHETA}.depth));

            observations.THETA{cTHETA}.depthIndex = ziTHETA';
            observations.THETA{cTHETA}.lon = (ones(1,length(THETA_obs)) .* rawObs{i}.lon(j))';
            observations.THETA{cTHETA}.lat = (ones(1,length(THETA_obs)) .* rawObs{i}.lat(j))';
            observations.THETA{cTHETA}.gridIndex = (ones(1,length(THETA_obs)) .* rawObs{i}.gridIndex(j))';
            observations.THETA{cTHETA}.type = (ones(1,length(THETA_obs)) .* rawObs{i}.type(j))';
            
            model.THETA{cTHETA}.data = THETA_model;
            
            cTHETA = cTHETA + 1;
            
        end
        
        %%
        %SAL
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziSAL,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziSAL(hi) = [];
        
        if (~isempty(ziSAL))
            
            SAL_obs = rawObs{i}.SAL(j,ziSAL)';
            
            SAL_path = [files.SAL(rawObs{i}.fileIndex).folder '/' files.SAL(rawObs{i}.fileIndex).name];
            SAL_model = read_ijk(SAL_path,rawObs{i}.gridIndex(j),1,ziSAL,nx*ny,1,nz) + 35;
            
            observations.SAL{cSAL}.data = SAL_obs;
            observations.SAL{cSAL}.time = (ones(1,length(SAL_obs)) .* rawObs{i}.time(j))';
            observations.SAL{cSAL}.depth = rawObs{i}.depth(j,ziSAL)';
            observations.SAL{cSAL}.depthIndex = ziSAL';
            observations.SAL{cSAL}.lon = (ones(1,length(SAL_obs)) .* rawObs{i}.lon(j))';
            observations.SAL{cSAL}.lat = (ones(1,length(SAL_obs)) .* rawObs{i}.lat(j))';
            observations.SAL{cSAL}.gridIndex = (ones(1,length(SAL_obs)) .* rawObs{i}.gridIndex(j))';
            observations.SAL{cSAL}.type = (ones(1,length(SAL_obs)) .* rawObs{i}.type(j))';
            
            model.SAL{cSAL}.data = SAL_model;
            
            cSAL = cSAL + 1;
            
        end
        
        %%
        %NO3
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziNO3,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziNO3(hi) = [];
        
        if (~isempty(ziNO3))
            
            NO3_obs = rawObs{i}.NO3(j,ziNO3)';
            
            NO3_path = [files.NO3(rawObs{i}.fileIndex).folder '/' files.NO3(rawObs{i}.fileIndex).name];
            NO3_model = read_ijk(NO3_path,rawObs{i}.gridIndex(j),1,ziNO3,nx*ny,1,nz);
            
            observations.NO3{cNO3}.data = NO3_obs;
            observations.NO3{cNO3}.time = (ones(1,length(NO3_obs)) .* rawObs{i}.time(j))';
            observations.NO3{cNO3}.depth = rawObs{i}.depth(j,ziNO3)';
            observations.NO3{cNO3}.depthIndex = ziNO3';
            observations.NO3{cNO3}.lon = (ones(1,length(NO3_obs)) .* rawObs{i}.lon(j))';
            observations.NO3{cNO3}.lat = (ones(1,length(NO3_obs)) .* rawObs{i}.lat(j))';
            observations.NO3{cNO3}.gridIndex = (ones(1,length(NO3_obs)) .* rawObs{i}.gridIndex(j))';
            observations.NO3{cNO3}.type = (ones(1,length(NO3_obs)) .* rawObs{i}.type(j))';
            
            model.NO3{cNO3}.data = NO3_model;
            
            cNO3 = cNO3 + 1;
            
        end
        
        %%
        %NO2
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziNO2,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziNO2(hi) = [];
        
        if (~isempty(ziNO2))
            
            NO2_obs = rawObs{i}.NO2(j,ziNO2)';
            
            NO2_path = [files.NO2(rawObs{i}.fileIndex).folder '/' files.NO2(rawObs{i}.fileIndex).name];
            NO2_model = read_ijk(NO2_path,rawObs{i}.gridIndex(j),1,ziNO2,nx*ny,1,nz);
            
            observations.NO2{cNO2}.data = NO2_obs;
            observations.NO2{cNO2}.time = (ones(1,length(NO2_obs)) .* rawObs{i}.time(j))';
            observations.NO2{cNO2}.depth = rawObs{i}.depth(j,ziNO2)';
            observations.NO2{cNO2}.depthIndex = ziNO2';
            observations.NO2{cNO2}.lon = (ones(1,length(NO2_obs)) .* rawObs{i}.lon(j))';
            observations.NO2{cNO2}.lat = (ones(1,length(NO2_obs)) .* rawObs{i}.lat(j))';
            observations.NO2{cNO2}.gridIndex = (ones(1,length(NO2_obs)) .* rawObs{i}.gridIndex(j))';
            observations.NO2{cNO2}.type = (ones(1,length(NO2_obs)) .* rawObs{i}.type(j))';
            
            model.NO2{cNO2}.data = NO2_model;
            
            cNO2 = cNO2 + 1;
            
        end
        
        %%
        %PO4
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziPO4,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziPO4(hi) = [];
        
        if (~isempty(ziPO4))
            
            PO4_obs = rawObs{i}.PO4(j,ziPO4)';
            
            PO4_path = [files.PO4(rawObs{i}.fileIndex).folder '/' files.PO4(rawObs{i}.fileIndex).name];
            PO4_model = read_ijk(PO4_path,rawObs{i}.gridIndex(j),1,ziPO4,nx*ny,1,nz);
            
            observations.PO4{cPO4}.data = PO4_obs;
            observations.PO4{cPO4}.time = (ones(1,length(PO4_obs)) .* rawObs{i}.time(j))';
            observations.PO4{cPO4}.depth = rawObs{i}.depth(j,ziPO4)';
            observations.PO4{cPO4}.depthIndex = ziPO4';
            observations.PO4{cPO4}.lon = (ones(1,length(PO4_obs)) .* rawObs{i}.lon(j))';
            observations.PO4{cPO4}.lat = (ones(1,length(PO4_obs)) .* rawObs{i}.lat(j))';
            observations.PO4{cPO4}.gridIndex = (ones(1,length(PO4_obs)) .* rawObs{i}.gridIndex(j))';
            observations.PO4{cPO4}.type = (ones(1,length(PO4_obs)) .* rawObs{i}.type(j))';
            
            model.PO4{cPO4}.data = PO4_model;
            
            cPO4 = cPO4 + 1;
            
        end
        
        %%
        %SiO2
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziSiO2,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziSiO2(hi) = [];
        
        if (~isempty(ziSiO2))
            
            SiO2_obs = rawObs{i}.SiO2(j,ziSiO2)';
            
            SiO2_path = [files.SiO2(rawObs{i}.fileIndex).folder '/' files.SiO2(rawObs{i}.fileIndex).name];
            SiO2_model = read_ijk(SiO2_path,rawObs{i}.gridIndex(j),1,ziSiO2,nx*ny,1,nz);
            
            observations.SiO2{cSiO2}.data = SiO2_obs;
            observations.SiO2{cSiO2}.time = (ones(1,length(SiO2_obs)) .* rawObs{i}.time(j))';
            observations.SiO2{cSiO2}.depth = rawObs{i}.depth(j,ziSiO2)';
            observations.SiO2{cSiO2}.depthIndex = ziSiO2';
            observations.SiO2{cSiO2}.lon = (ones(1,length(SiO2_obs)) .* rawObs{i}.lon(j))';
            observations.SiO2{cSiO2}.lat = (ones(1,length(SiO2_obs)) .* rawObs{i}.lat(j))';
            observations.SiO2{cSiO2}.gridIndex = (ones(1,length(SiO2_obs)) .* rawObs{i}.gridIndex(j))';
            observations.SiO2{cSiO2}.type = (ones(1,length(SiO2_obs)) .* rawObs{i}.type(j))';
            
            model.SiO2{cSiO2}.data = SiO2_model;
            
            cSiO2 = cSiO2 + 1;
            
        end
        
        %%
        %O2
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziO2,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziO2(hi) = [];
        
        if (~isempty(ziO2))
            
            O2_obs = rawObs{i}.O2(j,ziO2)';
            
            O2_path = [files.O2(rawObs{i}.fileIndex).folder '/' files.O2(rawObs{i}.fileIndex).name];
            O2_model = read_ijk(O2_path,rawObs{i}.gridIndex(j),1,ziO2,nx*ny,1,nz);
            
            observations.O2{cO2}.data = O2_obs;
            observations.O2{cO2}.time = (ones(1,length(O2_obs)) .* rawObs{i}.time(j))';
            observations.O2{cO2}.depth = rawObs{i}.depth(j,ziO2)';
            observations.O2{cO2}.depthIndex = ziO2';
            observations.O2{cO2}.lon = (ones(1,length(O2_obs)) .* rawObs{i}.lon(j))';
            observations.O2{cO2}.lat = (ones(1,length(O2_obs)) .* rawObs{i}.lat(j))';
            observations.O2{cO2}.gridIndex = (ones(1,length(O2_obs)) .* rawObs{i}.gridIndex(j))';
            observations.O2{cO2}.type = (ones(1,length(O2_obs)) .* rawObs{i}.type(j))';
            
            model.O2{cO2}.data = O2_model;
            
            cO2 = cO2 + 1;
            
        end
        
        %%
        %DOC
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziDOC,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziDOC(hi) = [];
        
        if (~isempty(ziDOC))
            
            DOC_obs = rawObs{i}.DOC(j,ziDOC)';
            
            DOC_path = [files.DOC(rawObs{i}.fileIndex).folder '/' files.DOC(rawObs{i}.fileIndex).name];
            DOC_model = read_ijk(DOC_path,rawObs{i}.gridIndex(j),1,ziDOC,nx*ny,1,nz);
            
            observations.DOC{cDOC}.data = DOC_obs;
            observations.DOC{cDOC}.time = (ones(1,length(DOC_obs)) .* rawObs{i}.time(j))';
            observations.DOC{cDOC}.depth = rawObs{i}.depth(j,ziDOC)';
            observations.DOC{cDOC}.depthIndex = ziDOC';
            observations.DOC{cDOC}.lon = (ones(1,length(DOC_obs)) .* rawObs{i}.lon(j))';
            observations.DOC{cDOC}.lat = (ones(1,length(DOC_obs)) .* rawObs{i}.lat(j))';
            observations.DOC{cDOC}.gridIndex = (ones(1,length(DOC_obs)) .* rawObs{i}.gridIndex(j))';
            observations.DOC{cDOC}.type = (ones(1,length(DOC_obs)) .* rawObs{i}.type(j))';
            
            model.DOC{cDOC}.data = DOC_model;
            
            cDOC = cDOC + 1;
            
        end
        
        %%
        %DON
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziDON,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziDON(hi) = [];
        
        if (~isempty(ziDON))
            
            DON_obs = rawObs{i}.DON(j,ziDON)';
            
            DON_path = [files.DON(rawObs{i}.fileIndex).folder '/' files.DON(rawObs{i}.fileIndex).name];
            DON_model = read_ijk(DON_path,rawObs{i}.gridIndex(j),1,ziDON,nx*ny,1,nz);
            
            observations.DON{cDON}.data = DON_obs;
            observations.DON{cDON}.time = (ones(1,length(DON_obs)) .* rawObs{i}.time(j))';
            observations.DON{cDON}.depth = rawObs{i}.depth(j,ziDON)';
            observations.DON{cDON}.depthIndex = ziDON';
            observations.DON{cDON}.lon = (ones(1,length(DON_obs)) .* rawObs{i}.lon(j))';
            observations.DON{cDON}.lat = (ones(1,length(DON_obs)) .* rawObs{i}.lat(j))';
            observations.DON{cDON}.gridIndex = (ones(1,length(DON_obs)) .* rawObs{i}.gridIndex(j))';
            observations.DON{cDON}.type = (ones(1,length(DON_obs)) .* rawObs{i}.type(j))';
            
            model.DON{cDON}.data = DON_model;
            
            cDON = cDON + 1;
            
        end
        
        %%
        %DIC
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziDIC,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziDIC(hi) = [];
        
        if (~isempty(ziDIC))
            
            DIC_obs = rawObs{i}.DIC(j,ziDIC)';
            DIC_pH_25C_obs = rawObs{i}.DIC_pH_25C(j,ziDIC)';
            DIC_pH_insitu_obs = rawObs{i}.DIC_pH_insitu(j,ziDIC)';
            
            DIC_path = [files.DIC(rawObs{i}.fileIndex).folder '/' files.DIC(rawObs{i}.fileIndex).name];
            DIC_model = read_ijk(DIC_path,rawObs{i}.gridIndex(j),1,ziDIC,nx*ny,1,nz);
            
            observations.DIC{cDIC}.data = DIC_obs;
            observations.DIC{cDIC}.data_pH_25C = DIC_pH_25C_obs;
            observations.DIC{cDIC}.data_pH_insitu = DIC_pH_insitu_obs;
            
            observations.DIC{cDIC}.time = (ones(1,length(DIC_obs)) .* rawObs{i}.time(j))';
            observations.DIC{cDIC}.depth = rawObs{i}.depth(j,ziDIC)';
            observations.DIC{cDIC}.depthIndex = ziDIC';
            observations.DIC{cDIC}.lon = (ones(1,length(DIC_obs)) .* rawObs{i}.lon(j))';
            observations.DIC{cDIC}.lat = (ones(1,length(DIC_obs)) .* rawObs{i}.lat(j))';
            observations.DIC{cDIC}.gridIndex = (ones(1,length(DIC_obs)) .* rawObs{i}.gridIndex(j))';
            observations.DIC{cDIC}.type = (ones(1,length(DIC_obs)) .* rawObs{i}.type(j))';
            
            model.DIC{cDIC}.data = DIC_model;
            
            cDIC = cDIC + 1;
            
        end
        
        %%
        %ALK
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziALK,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziALK(hi) = [];
        
        if (~isempty(ziALK))
            
            ALK_obs = rawObs{i}.ALK(j,ziALK)';
            
            ALK_path = [files.ALK(rawObs{i}.fileIndex).folder '/' files.ALK(rawObs{i}.fileIndex).name];
            ALK_model = read_ijk(ALK_path,rawObs{i}.gridIndex(j),1,ziALK,nx*ny,1,nz);
            
            observations.ALK{cALK}.data = ALK_obs;
            observations.ALK{cALK}.time = (ones(1,length(ALK_obs)) .* rawObs{i}.time(j))';
            observations.ALK{cALK}.depth = rawObs{i}.depth(j,ziALK)';
            observations.ALK{cALK}.depthIndex = ziALK';
            observations.ALK{cALK}.lon = (ones(1,length(ALK_obs)) .* rawObs{i}.lon(j))';
            observations.ALK{cALK}.lat = (ones(1,length(ALK_obs)) .* rawObs{i}.lat(j))';
            observations.ALK{cALK}.gridIndex = (ones(1,length(ALK_obs)) .* rawObs{i}.gridIndex(j))';
            observations.ALK{cALK}.type = (ones(1,length(ALK_obs)) .* rawObs{i}.type(j))';
            
            model.ALK{cALK}.data = ALK_model;
            
            cALK = cALK + 1;
            
        end
        
        %%
        %PH
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziPH,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziPH(hi) = [];
        
        if (~isempty(ziPH))
            
            PH_obs = rawObs{i}.PH(j,ziPH)';
            
            PH_path = [files.PH(rawObs{i}.fileIndex).folder '/' files.PH(rawObs{i}.fileIndex).name];
            PH_model = read_ijk(PH_path,rawObs{i}.gridIndex(j),1,ziPH,nx*ny,1,nz);
            
            observations.PH{cPH}.data = PH_obs;
            observations.PH{cPH}.time = (ones(1,length(PH_obs)) .* rawObs{i}.time(j))';
            observations.PH{cPH}.depth = rawObs{i}.depth(j,ziPH)';
            observations.PH{cPH}.depthIndex = ziPH';
            observations.PH{cPH}.lon = (ones(1,length(PH_obs)) .* rawObs{i}.lon(j))';
            observations.PH{cPH}.lat = (ones(1,length(PH_obs)) .* rawObs{i}.lat(j))';
            observations.PH{cPH}.gridIndex = (ones(1,length(PH_obs)) .* rawObs{i}.gridIndex(j))';
            observations.PH{cPH}.type = (ones(1,length(PH_obs)) .* rawObs{i}.type(j))';
            
            model.PH{cPH}.data = PH_model;
            
            cPH = cPH + 1;
            
        end
        
        %%
        %CHLA
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
        
        hi = find(hFacC == 0);
        
        ziCHLA(hi) = [];
        
        if (~isempty(ziCHLA))
            
            CHLA_obs = rawObs{i}.CHLA(j,ziCHLA)';
            
            CHL_path1 = [files.CHL1(rawObs{i}.fileIndex).folder '/' files.CHL1(rawObs{i}.fileIndex).name];
            CHL_path2 = [files.CHL2(rawObs{i}.fileIndex).folder '/' files.CHL2(rawObs{i}.fileIndex).name];
            CHL_path3 = [files.CHL3(rawObs{i}.fileIndex).folder '/' files.CHL3(rawObs{i}.fileIndex).name];
            CHL_path4 = [files.CHL4(rawObs{i}.fileIndex).folder '/' files.CHL4(rawObs{i}.fileIndex).name];
            CHL_path5 = [files.CHL5(rawObs{i}.fileIndex).folder '/' files.CHL5(rawObs{i}.fileIndex).name];
            
            CHL1_model = read_ijk(CHL_path1,rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
            CHL2_model = read_ijk(CHL_path2,rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
            CHL3_model = read_ijk(CHL_path3,rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
            CHL4_model = read_ijk(CHL_path4,rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
            CHL5_model = read_ijk(CHL_path5,rawObs{i}.gridIndex(j),1,ziCHLA,nx*ny,1,nz);
            
            CHL1_model(CHL1_model < 0) = 0;
            CHL2_model(CHL2_model < 0) = 0;
            CHL3_model(CHL3_model < 0) = 0;
            CHL4_model(CHL4_model < 0) = 0;
            CHL5_model(CHL5_model < 0) = 0;
            
            CHL_model = CHL1_model + CHL2_model + CHL3_model + CHL4_model + CHL5_model;
            
            observations.CHLA{cCHLA}.data = CHLA_obs;
            observations.CHLA{cCHLA}.time = (ones(1,length(CHLA_obs)) .* rawObs{i}.time(j))';
            observations.CHLA{cCHLA}.depth = rawObs{i}.depth(j,ziCHLA)';
            observations.CHLA{cCHLA}.depthIndex = ziCHLA';
            observations.CHLA{cCHLA}.lon = (ones(1,length(CHLA_obs)) .* rawObs{i}.lon(j))';
            observations.CHLA{cCHLA}.lat = (ones(1,length(CHLA_obs)) .* rawObs{i}.lat(j))';
            observations.CHLA{cCHLA}.gridIndex = (ones(1,length(CHLA_obs)) .* rawObs{i}.gridIndex(j))';
            observations.CHLA{cCHLA}.type = (ones(1,length(CHLA_obs)) .* rawObs{i}.type(j))';
            
            model.CHLA{cCHLA}.data = CHL_model;
            
            cCHLA = cCHLA + 1;
            
        end
        
    end
    
    disp(num2str(i));
    fprintf(fileID,'file: %d \n',i);
    
end

%%

save([saveDir 'GLODAP_model_data_pairs.mat'],'model','observations','-v7.3');

fclose(fileID);

cd(codeDir);

end

%%
