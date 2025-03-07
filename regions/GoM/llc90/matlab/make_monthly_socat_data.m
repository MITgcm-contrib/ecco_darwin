
%This code takes the downloaded socat fCO2 data and turns it into daily
%means as most days have multiple entries. 

addpath('C:\Users\jessicaz\Documents\MATLAB\lib')

load C:\Users\jessicaz\Documents\research\socat\gomscat.mat


%{ 
% for kk=1:length(yr)
% 
%     dates(kk)=datetime([yr{kk} mon{kk} day{kk}],'InputFormat','yyyyMMdd');
% 
% end
% 
% save('C:\Users\jessicaz\Documents\research\socat\gomscat.mat','dates','-append')
% 

% latitude=cellfun(@str2num,latitude(:),'un',0);
% longitude=cellfun(@str2num,longitude(:),'un',0);
% 
% latitude=[latitude{:}]';
% longitude=[longitude{:}]';
% 
% save('C:\Users\jessicaz\Documents\research\socat\gomscat.mat',...
%     'longitude','latitude','-append')
%}

nx=20;
ny=15;

%Read in the model grid to use to discretize the socat data
mlat=readbin('C:\Users\jessicaz\Documents\research\regions\GoM\llc90\grid\YG_20x15',...
    [nx ny]);
mlat=mlat(1,:)';

mlon=readbin('C:\Users\jessicaz\Documents\research\regions\GoM\llc90\grid\XG_20x15',...
    [nx ny]);
mlon=mlon(:,1);


%Distrectize socat data into llc90 grid
[N,~,~,binX,binY] = histcounts2(longitude,latitude,...
    [mlon; mlon(end)+mlon(end)-mlon(end-1)],...
    [mlat; mlat(end)+mlat(end)-mlat(end-1)]);

%Create array with the xbin, ybin, date, and mean fCO2 for that date for
%the socat data

R=unique(dates.Year); 

fco2rec=cellfun(@str2num,fco2rec(:),'un',0);
fco2rec=[fco2rec{:}]';

sal=cellfun(@str2num,sal(:),'un',0);
sal=[sal{:}]';

sst=cellfun(@str2num,sst(:),'un',0);
sst=[sst{:}]';

%make data storage arrays
gom_fco2=nan(nx,ny);
gom_sst=nan(nx,ny);
gom_sal=nan(nx,ny);

allDate=datetime('today');
for rr=1:length(R) %loop through years

    for tmpMon=1:12 %loop through months
    
        datatmp1=nan(nx,ny);
        datatmp2=nan(nx,ny);
        datatmp3=nan(nx,ny);

        dind=find(dates.Year == R(rr) & dates.Month == tmpMon);

        tmpX=binX(dind);
        tmpY=binY(dind);
        tmpfCO2=fco2rec(dind);
        tmpsst=sst(dind);
        tmpsal=sal(dind);

       %loop through xy grid
      for xx=1:nx
            for yy=1:ny
    
                ind=find(tmpX==xx & tmpY==yy);
    
                datatmp1(xx,yy)=mean(tmpfCO2(ind));
                datatmp2(xx,yy)=mean(tmpsst(ind));
                datatmp3(xx,yy)=mean(tmpsal(ind));    
    
            end
      end
    
        gom_fco2=cat(3,gom_fco2,datatmp1);
        gom_sst=cat(3,gom_sst,datatmp2);
        gom_sal=cat(3,gom_sal,datatmp3);

        allDate=cat(1,allDate,datetime([R(rr) tmpMon 1]));

    end
    
end

save('gom_grid_data.mat','gom_sal','gom_sst','gom_fco2','mlat','mlon',...
    'R','allDate')

figure()
kk=8;

pcolorcen(gom_fco2(:,:,kk));
colorbar()
kk=kk+1
