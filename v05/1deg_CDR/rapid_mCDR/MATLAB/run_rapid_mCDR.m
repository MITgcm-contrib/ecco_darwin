clear
close all;

%%

dataDir = '/Users/carrolld/Documents/research/mCDR/rapid_mCDR/input/';

expName = {'NAS', 'WBC', 'EBC', 'EU', 'ACC'};

avgType = ''; %set to '_pco2' for HorAdv averaging, leave as '' for Deploy averaging

%%

%function to compute cumulative mean
cMean = @(dta) cumsum(dta, 'omitnan') ./ (1:numel(dta));

out = struct();

%%

for i = 1:length(expName)
    
    exp = expName{i};
    
    path = [dataDir, exp, '/model_forcing_data_daily', avgType, '.nc'];
    
    ni = ncinfo(path);
    
    %ncdisp(path);
    
    for j = 1:length(ni.Variables)
        
        vn = ni.Variables(j).Name;
        input.(vn) = ncread(path, vn)';
        
    end

    input.area2 = repmat(input.area,[(length(input.Zl)-1) 1])';
    
    out.(exp) = rapid_mCDR_model(input);

    %%
    
    out.(exp).time = double(out.(exp).time/86400) + datenum(1992,1,1,0,0,0);

    out.(exp).X = (-out.(exp).fCO2 ./ out.(exp).alk_forcing' .* out.(exp).area') ./ ...
        cMean(out.(exp).cdr_potential)';
    
    [years months days hours minutes seconds] = datevec(out.(exp).time);
    
    uniqueYears = unique(years);
    
    for k = 1:length(uniqueYears)
        
        ti = find(years == uniqueYears(k));
        
        out.(exp).yearlyTime(k) = datenum(uniqueYears(k),6,15,0,0,0);

        out.(exp).yearlyFCO2(k) = nanmean(out.(exp).fCO2(ti));
        out.(exp).yearlyX(k) = nanmean(out.(exp).X(ti));
        
    end
    
    disp(exp);
    
end

%%

years = 1995:2017;

dateTicks = [datenum(1995,1,1,0,0,0) datenum(2000,1,1,0,0,0) datenum(2005,1,1,0,0,0) ...
    datenum(2010,1,1,0,0,0) datenum(2015,1,1,0,0,0) datenum(2020,1,1,0,0,0)];

startDate = datenum(years(1),1,1,0,0,0);
endDate = datenum(years(end),12,31,0,0,0);

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

lw = 2;
fs = 26;

wind = 365;

colors1 = cbrewer('qual','Set1',6);

colors2(1,:) = colors1(2,:);
colors2(2,:) = colors1(5,:);
colors2(3,:) = colors1(3,:);
colors2(4,:) = colors1(1,:);
colors2(5,:) = colors1(4,:);

subplot(211);

hold on

for i = 1:length(expName)
    
    exp = expName{i};
    
    eval(['p' num2str(i) '= plot(out.(exp).time,movmean(out.(exp).fCO2, wind,''Endpoints'',''fill''),''Color'',colors2(i,:),''LineWidth'',1);']);
    eval(['p' num2str(i) '= plot(out.(exp).yearlyTime,out.(exp).yearlyFCO2,''Color'',colors2(i,:),''LineWidth'',lw);']);

end

legend([p1 p2 p3 p4 p5],expName,'Orientation','Horizontal','Location','Best');

axis tight

ylabel('CO_2 Flux (mmol C m^-^2)');

set(gca,'xtick',dateTicks);

datetick('x',11,'keepticks');

xlim([startDate endDate]);

box on
grid on
set(gca,'GridLineStyle','-.');
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title('CO_2 Flux');

subplot(212);

hold on

for i = 1:length(expName)
    
    exp = expName{i};
    
    eval(['p' num2str(i) '=  plot(out.(exp).time,movmean(out.(exp).X, wind,''Endpoints'',''fill''),''Color'',colors2(i,:),''LineWidth'',lw);']);
    eval(['p' num2str(i) '= plot(out.(exp).yearlyTime,out.(exp).yearlyX,''Color'',colors2(i,:),''LineWidth'',lw);']);

end

xlabel('Year');
ylabel('CDR Efficiency');

set(gca,'xtick',dateTicks);
set(gca,'ytick',[0:0.2:2]);

datetick('x',11,'keepticks');

xlim([startDate endDate]);
ylim([0 1.4]);

box on
grid on
set(gca,'GridLineStyle','-.');
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title('CDR Efficiency');

%%
