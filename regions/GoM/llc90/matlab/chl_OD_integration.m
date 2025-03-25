%This code makes binary files of the chlorophyll integrated over the
%optical depth range as defined by globcolor data

addpath('C:\Users\jessicaz\Documents\MATLAB\lib')

load GlobColor\data\kd490_20x15x328.mat
load GlobColor\data\chl_20x15x328.mat


%for the months with globcolor data, find the corresponding llc90 output.
%then integrate over the depth as defined by the kd490. then compare to the
%globcolor chl data that is also regridded







