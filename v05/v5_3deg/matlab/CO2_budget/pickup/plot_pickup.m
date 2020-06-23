clear
close all

dataDir = '/Users/carrolld/Documents/research/v4_3deg_JAMES_budget_latest/pickup/';

field = readbin([dataDir 'ptracers_optimized_01.0000000001'],[64 128 15],1,'real*8');

pcolorcen(field(:,:,1));