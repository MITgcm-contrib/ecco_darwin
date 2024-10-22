clear
close all;

%%

load(['bin_average_LLC_270_to_1x1_deg.mat']);

%%

%LLC 270 grid dimensions
nx1 = 270;
ny1 = nx1 .* 13;

%1x1 deg bin-averaged grid dimensions
nx2 = 360;
ny2 = 180;

%read raw CO2 flux field
CO2_flux = -readbin(['DICCFLXmonth.0000002232.data'],[nx1 ny1]);

%%
%bin-averaged grid

lon = output.XG;
lat = output.YG;
area = output.RAC;

%%

temp1 = input.hFacC;
temp2 = CO2_flux;

output.hFacC = reshape(bin_average_cons*temp1(:),[nx2 ny2]);

%bin average using flux conservation
output.binField = reshape(bin_average_cons*temp2(:),[nx2 ny2]);

%if computing non-flux fields use the following
%output.binField = reshape(bin_average_cons*temp2(:),[nx2 ny2]) ./ output.hFacC;

%mask dry cells for plotting purposes
output.binField(output.hFacC == 0) = nan;

%%

hFig1 = figure(1);
set(hFig1,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color',[1 1 1]);

lw = 2;
fs = 26;

cMin = -5.*10^-4;
cMax = 5.*10^-4;

hold on

set(gca,'Color',[0.6 0.6 0.6]);

pcolorcen(lon,lat,output.binField);

caxis([cMin cMax]);

hcb = colorbar;

set(get(hcb,'ylabel'),'String',{'Air-sea CO_2 Flux (mmol m^-^2 s^-^1)'}');

axis tight

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

box on
set(gca,'LineWidth',lw);
set(gca,'FontSize',fs);

title('Bin-averaged Air-sea CO_2 Flux');

%%
%test flux conservation

intCO2FluxRaw = nansum(CO2_flux(:) .* input.RAC(:)); %globally-integrated CO2 flux from raw field
intCO2FluxBin = nansum(output.binField(:) .* area(:)); %globally-integrated CO2 flux from bin-averaged field

%%
