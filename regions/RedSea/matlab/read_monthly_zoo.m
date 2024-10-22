clear all

%%%% This routine reads monthly fields of zooplankton %%%%%%%%%%%%%
%%%% It should be placed in the 'run' directory of the run that will be processed
%%%% cp /ecco_darwin/regions/RedSea/matlab/read_monthly_zoo.m /darwin3/run/diags/monthly/ 


zoo_micro_monthly = rdmds('c6',NaN); %%%%% read all months of directory (mmol C   m^-3)
zoo_micro_monthly_sur=squeeze(zoo_micro_monthly(:,:,1,:));%%%%% pick depth, 1 for surface
zoo_micro_monthly_sur_ave=squeeze(nanmean(nanmean(zoo_micro_monthly_sur)));%%%%% average the whole basin 


%%%%% plot all months basin average
figure (1)
plot(zoo_micro_monthly_sur_ave,'linewidth',2)
ylabel('Zooplankton Concentration (mmol C   m^-^3)')
xlabel('Time (month)') %%%% first month is January 1992 %%%%%
set(gca,'fontsize',14)


%%%%% read model grid coordinates and depth %%%%%%
G=load_grid('../..');
X=G.xC;
Y=G.yC;


%%%%% Map a month's average 
figure(2)
pcolorcen(X,Y,zoo_micro_monthly_sur(:,:,1))




