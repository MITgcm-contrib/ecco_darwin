function [] = model_obs_compare_profiles(loadData,modelDir)

gridDir = '/nobackup/dcarrol2/grid/LLC_270/';

dataDir = '/nobackup/dcarrol2/evaluation/mat/observations/BGC-Argo/';
codeDir = '/nobackup/dcarrol2/evaluation/m_files/model_obs_compare/BGC-Argo/online/';
saveDir = '/nobackup/dcarrol2/evaluation/mat/model_obs_compare/BGC-Argo/';

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

load([dataDir 'BGC-Argo_profiles_ver_1.2.mat']);

fileList = dir([modelDir 'O2.0*.*data']);

for i = 1:length(fileList)
    
    fileName = fileList(i).name;
    
    temp = datenum(ts2dte(str2num(fileName(end-14:end-5)),1200,1992));

    modelTimeLowerBound(i) = addtodate(double(temp), -1, 'month');
    
    temp = addtodate(double(temp), -1, 'day');
    
    [year month day hour minute seconds] = datevec(temp);
    
    modelTimeUpperBound(i) = datenum(year,month,eomday(year,month),23,59,59);
    
end

%%

c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;

for i = 1:length(modelTimeLowerBound) %for all model output
    
    it1 = find(profiles.NO3.time >= modelTimeLowerBound(i) & profiles.NO3.time <= modelTimeUpperBound(i));
    it2 = find(profiles.O2.time >= modelTimeLowerBound(i) & profiles.O2.time <= modelTimeUpperBound(i));
    it3 = find(profiles.CHL.time >= modelTimeLowerBound(i) & profiles.CHL.time <= modelTimeUpperBound(i));
    it4 = find(profiles.PH.time >= modelTimeLowerBound(i) & profiles.PH.time <= modelTimeUpperBound(i));
    
    %%
    
    if (~isempty(it1))
        
        lon = profiles.NO3.lon(it1);
        lat = profiles.NO3.lat(it1);
        
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
        
        ix = find(landMask(im) == 1);
        
        rawObs1{c1}.time = profiles.NO3.time(it1(ix));
        rawObs1{c1}.lon = profiles.NO3.lon(it1(ix));
        rawObs1{c1}.lat = profiles.NO3.lat(it1(ix));
        rawObs1{c1}.type = profiles.NO3.type(it1(ix));
        
        rawObs1{c1}.numObs = length(im(ix));
        rawObs1{c1}.gridIndex = im(ix);
        rawObs1{c1}.fileIndex = i;
        
        rawObs1{c1}.NO3 = profiles.NO3.NO3(it1(ix),:);
        rawObs1{c1}.depth = profiles.NO3.depth(it1(ix),:);

        c1 = c1 + 1;
    
	clear ix 
     
    end
    
    clear im
    
    %%

    if (~isempty(it2))
        
        lon = profiles.O2.lon(it2);
        lat = profiles.O2.lat(it2);
        
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
        
        ix = find(landMask(im) == 1);
        
        rawObs2{c2}.time = profiles.O2.time(it2(ix));
        rawObs2{c2}.lon = profiles.O2.lon(it2(ix));
        rawObs2{c2}.lat = profiles.O2.lat(it2(ix));
        rawObs2{c2}.type = profiles.O2.type(it2(ix));
        
        rawObs2{c2}.numObs = length(im(ix));
        rawObs2{c2}.gridIndex = im(ix);
        rawObs2{c2}.fileIndex = i;
        
        rawObs2{c2}.O2 = profiles.O2.O2(it2(ix),:);
        rawObs2{c2}.depth = profiles.O2.depth(it2(ix),:); 
        
	c2 = c2 + 1;

	clear ix
        
    end
    
    clear im
    
    %%
    
    if (~isempty(it3))
        
        lon = profiles.CHL.lon(it3);
        lat = profiles.CHL.lat(it3);
        
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
        
        ix = find(landMask(im) == 1);
        
        rawObs3{c3}.time = profiles.CHL.time(it3(ix));
        rawObs3{c3}.lon = profiles.CHL.lon(it3(ix));
        rawObs3{c3}.lat = profiles.CHL.lat(it3(ix));
        rawObs3{c3}.type = profiles.CHL.type(it3(ix));
        
        rawObs3{c3}.numObs = length(im(ix));
        rawObs3{c3}.gridIndex = im(ix);
        rawObs3{c3}.fileIndex = i;
        
        rawObs3{c3}.CHL = profiles.CHL.CHL(it3(ix),:);
        rawObs3{c3}.depth = profiles.CHL.depth(it3(ix),:);

	clear ix
        
        c3 = c3 + 1;
        
    end
    
    clear im
    
    %%

    if (~isempty(it4))
        
        lon = profiles.PH.lon(it4);
        lat = profiles.PH.lat(it4);
        
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
        
        ix = find(landMask(im) == 1);
        
        rawObs4{c4}.time = profiles.PH.time(it4(ix));
        rawObs4{c4}.lon = profiles.PH.lon(it4(ix));
        rawObs4{c4}.lat = profiles.PH.lat(it4(ix));
        rawObs4{c4}.type = profiles.PH.type(it4(ix));
        
        rawObs4{c4}.numObs = length(im(ix));
        rawObs4{c4}.gridIndex = im(ix);
        rawObs4{c4}.fileIndex = i;
        
        rawObs4{c4}.PH = profiles.PH.PH(it4(ix),:);
        rawObs4{c4}.depth = profiles.PH.depth(it4(ix),:);

        c4 = c4 + 1;
    
	clear ix    
        
    end
    
    clear im

    disp(num2str(i));
    
end

%%
%model output

files.NO3 = dir([modelDir 'NO3.0*.*data']);
files.O2 = dir([modelDir 'O2.0*.*data']);
files.PH = dir([modelDir 'pH.0*.*data']);
files.CHL1 = dir([modelDir 'Chl1.0*.*data']);
files.CHL2 = dir([modelDir 'Chl2.0*.*data']);
files.CHL3 = dir([modelDir 'Chl3.0*.*data']);
files.CHL4 = dir([modelDir 'Chl4.0*.*data']);
files.CHL5 = dir([modelDir 'Chl5.0*.*data']);

%%

cNO3 = 1;


disp(['NO3: ' num2str(length(rawObs1))]);
disp(['O2: ' num2str(length(rawObs2))]);
disp(['CHL: ' num2str(length(rawObs3))]);
disp(['PH: ' num2str(length(rawObs4))]);

for i = 1:length(rawObs1)
    
    for j = 1:rawObs1{i}.numObs 

        ziNO3 = find(~isnan(rawObs1{i}.NO3(j,:)));
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs1{i}.gridIndex(j),1,ziNO3,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziNO3(hi) = [];
        
        if (~isempty(ziNO3))
            
            NO3_obs = rawObs1{i}.NO3(j,ziNO3)';
            
            NO3_path = [files.NO3(rawObs1{i}.fileIndex).folder '/' files.NO3(rawObs1{i}.fileIndex).name];
            NO3_model = read_ijk(NO3_path,rawObs1{i}.gridIndex(j),1,ziNO3,nx*ny,1,nz);
            
            observations.NO3{cNO3}.data = NO3_obs;
            observations.NO3{cNO3}.time = (ones(1,length(NO3_obs)) .* rawObs1{i}.time(j))';
            observations.NO3{cNO3}.depth = rawObs1{i}.depth(j,ziNO3)';
    	    observations.NO3{cNO3}.depthIndex = ziNO3';
            observations.NO3{cNO3}.lon = (ones(1,length(NO3_obs)) .* rawObs1{i}.lon(j))';
            observations.NO3{cNO3}.lat = (ones(1,length(NO3_obs)) .* rawObs1{i}.lat(j))';
            observations.NO3{cNO3}.gridIndex = (ones(1,length(NO3_obs)) .* rawObs1{i}.gridIndex(j))';
            observations.NO3{cNO3}.type = (ones(1,length(NO3_obs)) .* rawObs1{i}.type(j))';
            
            model.NO3{cNO3}.data = NO3_model;
            
            cNO3 = cNO3 + 1;
            
        end
        
	clear ziNO3

    end

    disp(num2str(i));
    
end

%% 

cO2 = 1;

for i = 1:length(rawObs2)
    
    for j = 1:rawObs2{i}.numObs 

        ziO2 = find(~isnan(rawObs2{i}.O2(j,:)));
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs2{i}.gridIndex(j),1,ziO2,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziO2(hi) = [];
        
        if (~isempty(ziO2))
            
            O2_obs = rawObs2{i}.O2(j,ziO2)';
            
            O2_path = [files.O2(rawObs2{i}.fileIndex).folder '/' files.O2(rawObs2{i}.fileIndex).name];
            O2_model = read_ijk(O2_path,rawObs2{i}.gridIndex(j),1,ziO2,nx*ny,1,nz);
            
            observations.O2{cO2}.data = O2_obs;
            observations.O2{cO2}.time = (ones(1,length(O2_obs)) .* rawObs2{i}.time(j))';
            observations.O2{cO2}.depth = rawObs2{i}.depth(j,ziO2)';
    	    observations.O2{cO2}.depthIndex = ziO2';
            observations.O2{cO2}.lon = (ones(1,length(O2_obs)) .* rawObs2{i}.lon(j))';
            observations.O2{cO2}.lat = (ones(1,length(O2_obs)) .* rawObs2{i}.lat(j))';
            observations.O2{cO2}.gridIndex = (ones(1,length(O2_obs)) .* rawObs2{i}.gridIndex(j))';
            observations.O2{cO2}.type = (ones(1,length(O2_obs)) .* rawObs2{i}.type(j))';
            
            model.O2{cO2}.data = O2_model;
            
            cO2 = cO2 + 1;
            
        end
        
	clear ziO2

    end
   
    disp(num2str(i));

end

%% 

cCHL = 1;

for i = 1:length(rawObs3)
    
    for j = 1:rawObs3{i}.numObs 

        ziCHL = find(~isnan(rawObs3{i}.CHL(j,:)));
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziCHL(hi) = [];
        
        if (~isempty(ziCHL))
            
            CHL_obs = rawObs3{i}.CHL(j,ziCHL)';
      
            CHL_path1 = [files.CHL1(rawObs3{i}.fileIndex).folder '/' files.CHL1(rawObs3{i}.fileIndex).name];
            CHL_path2 = [files.CHL2(rawObs3{i}.fileIndex).folder '/' files.CHL2(rawObs3{i}.fileIndex).name];
            CHL_path3 = [files.CHL3(rawObs3{i}.fileIndex).folder '/' files.CHL3(rawObs3{i}.fileIndex).name];
            CHL_path4 = [files.CHL4(rawObs3{i}.fileIndex).folder '/' files.CHL4(rawObs3{i}.fileIndex).name];
            CHL_path5 = [files.CHL5(rawObs3{i}.fileIndex).folder '/' files.CHL5(rawObs3{i}.fileIndex).name];
            
            CHL1_model = read_ijk(CHL_path1,rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
            CHL2_model = read_ijk(CHL_path2,rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
            CHL3_model = read_ijk(CHL_path3,rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
            CHL4_model = read_ijk(CHL_path4,rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
            CHL5_model = read_ijk(CHL_path5,rawObs3{i}.gridIndex(j),1,ziCHL,nx*ny,1,nz);
            
            CHL_model = CHL1_model + CHL2_model + CHL3_model + CHL4_model + CHL5_model;
       
            observations.CHL{cCHL}.data = CHL_obs;
            observations.CHL{cCHL}.time = (ones(1,length(CHL_obs)) .* rawObs3{i}.time(j))';
            observations.CHL{cCHL}.depth = rawObs3{i}.depth(j,ziCHL)';
    	    observations.CHL{cCHL}.depthIndex = ziCHL';
            observations.CHL{cCHL}.lon = (ones(1,length(CHL_obs)) .* rawObs3{i}.lon(j))';
            observations.CHL{cCHL}.lat = (ones(1,length(CHL_obs)) .* rawObs3{i}.lat(j))';
            observations.CHL{cCHL}.gridIndex = (ones(1,length(CHL_obs)) .* rawObs3{i}.gridIndex(j))';
            observations.CHL{cCHL}.type = (ones(1,length(CHL_obs)) .* rawObs3{i}.type(j))';
            
            model.CHL{cCHL}.data = CHL_model;
            
            cCHL = cCHL + 1;
            
        end
        
	clear ziCHL

    end
    
    disp(num2str(i));

end

%% 

cPH = 1;

for i = 1:length(rawObs4)
    
    for j = 1:rawObs4{i}.numObs 

        ziPH = find(~isnan(rawObs4{i}.PH(j,:)));
        
        hFacC = read_ijk([gridDir 'hFacC.data'],rawObs4{i}.gridIndex(j),1,ziPH,nx*ny,1,nz);
        hi = find(hFacC == 0);
        
        ziPH(hi) = [];
        
        if (~isempty(ziPH))
            
            PH_obs = rawObs4{i}.PH(j,ziPH)';
            
            PH_path = [files.PH(rawObs4{i}.fileIndex).folder '/' files.PH(rawObs4{i}.fileIndex).name];
            PH_model = read_ijk(PH_path,rawObs4{i}.gridIndex(j),1,ziPH,nx*ny,1,nz);
            
            observations.PH{cPH}.data = PH_obs;
            observations.PH{cPH}.time = (ones(1,length(PH_obs)) .* rawObs4{i}.time(j))';
            observations.PH{cPH}.depth = rawObs4{i}.depth(j,ziPH)';
    	    observations.PH{cPH}.depthIndex = ziPH';
            observations.PH{cPH}.lon = (ones(1,length(PH_obs)) .* rawObs4{i}.lon(j))';
            observations.PH{cPH}.lat = (ones(1,length(PH_obs)) .* rawObs4{i}.lat(j))';
            observations.PH{cPH}.gridIndex = (ones(1,length(PH_obs)) .* rawObs4{i}.gridIndex(j))';
            observations.PH{cPH}.type = (ones(1,length(PH_obs)) .* rawObs4{i}.type(j))';
            
            model.PH{cPH}.data = PH_model;
            
            cPH = cPH + 1;
            
        end
 	
	clear ziPH

    end

    disp(num2str(i));
    
end

%%

save([saveDir 'BGC-Argo_model_data_pairs.mat'],'model','observations','-v7.3');

%%
