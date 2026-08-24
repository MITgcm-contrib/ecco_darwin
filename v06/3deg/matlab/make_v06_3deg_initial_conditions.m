%% make_v06_3deg_initial_conditions.m
%
% Build the 36 ptracer initial-condition files needed by v06/3deg from the 31
% files that ship with v05/3deg/data_darwin.
%
% Why this is not a straight copy
% ------------------------------
% v05 and v06 have different ecosystems, so plankton tracer INDICES do not
% correspond to the same organisms:
%
%   v05 (nplank=7,  nPhoto=5)   c01 diatom (hasSi)
%                               c02 cocco  (hasPIC)
%                               c03 cocco  (hasPIC)
%                               c04 other phyto
%                               c05 other phyto
%                               c06 zoo
%                               c07 zoo
%
%   v06 (nplank=10, nPhoto=6)   c01 PicoCyano
%                               c02 PicoEuk
%                               c03 Cocco
%                               c04 Diazo      (absent from v05/3deg: DIAZO=7*0)
%                               c05 Diatom
%                               c06 Dino
%                               c07..c10 Zoo (4 size classes)
%
% There is no honest one-to-one map. What this script does instead is conserve
% total carbon: all v05 phytoplankton biomass is pooled and split evenly across
% the 6 v06 phytoplankton, all v05 zooplankton biomass is pooled and split
% evenly across the 4 v06 zooplankton. Darwin sorts out community structure
% within a few model months, so the split only sets the spin-up transient, not
% the equilibrium. Total plankton C is preserved exactly.
%
% Tracers 1-19 (DIC..O2) are copied unchanged -- these are the ones that
% actually matter for the initial state and they map 1:1 between versions.
% CDOM (20) starts at zero and fills in from DOC production (fracCDOM=0.02).
% Chl01..Chl06 (31-36) are NOT written: darwin_chlInitBalanced=T in data.darwin
% computes chlorophyll from biomass at darwin_chlIter0.
%
% A note on negatives
% -------------------
% The v05/3deg optimized initial conditions contain small negative values in
% most tracers (e.g. NO3 min = -3.13, O2 min = -27.2) -- an artifact of the
% Green's Functions optimization, present in v05 as shipped. v06 turns on
% DARWIN_ALLOW_CONS, so these will show up in the conservation monitor.
% Set clipNegatives=true below to floor them at zero; leave it false to
% reproduce the v05 initial state bit-for-bit in tracers 1-19.
%
% Usage: edit inDir/outDir below, then run.

%% paths and grid

inDir  = '/Users/carrolld/Documents/research/debug/ecco_darwin/v05/3deg/data_darwin/';
outDir = '/Users/carrolld/Documents/research/debug/v06/3deg/data_darwin/';

nx = 128;
ny =  64;
nz =  15;   % Nr in code/SIZE.h

% readBinaryPrec=32 in input/data -> real*4
fmt  = 'ieee-be';  % MITgcm default byte order

clipNegatives = false;   % see "A note on negatives" above

if ~exist(outDir,'dir'); mkdir(outDir); end

inSuff  = '.0000000001';
outSuff = '.0000000001';

%% 1-19: DIC, NO3, NO2, NH4, PO4, FeT, SiO2, DOC, DON, DOP, DOFe,
%%       POC, PON, POP, POFe, POSi, PIC, ALK, O2  -- direct copy

for n = 1:19
    src = fullfile(inDir, sprintf('ptracers_optimized_%02d%s', n, inSuff));
    dst = fullfile(outDir, sprintf('ptracers_v06_3deg_%02d%s', n, outSuff));

    fld = read3d(src, nx, ny, nz, fmt);
    if clipNegatives
        fld = max(fld, 0);
    end
    write3d(dst, fld, fmt);

    fprintf('%2d  %-6s  copied   min=%10.4g  max=%10.4g\n', ...
            n, tracerName(n), min(fld(:)), max(fld(:)));
end

%% 20: CDOM -- new in v06, start from zero

cdom = zeros(nx, ny, nz);
write3d(fullfile(outDir, sprintf('ptracers_v06_3deg_20%s', outSuff)), cdom, fmt);
fprintf('20  CDOM    zeros (new tracer, spins up from fracCDOM)\n');

%% 21-30: plankton -- pool v05 biomass, redistribute over the v06 groups

% v05 phytoplankton are tracers 20-24 (c01..c05), zooplankton are 25-26 (c06..c07)
v05phyto = 20:24;
v05zoo   = 25:26;

phytoSum = zeros(nx, ny, nz);
for n = v05phyto
    phytoSum = phytoSum + read3d(fullfile(inDir, ...
        sprintf('ptracers_optimized_%02d%s', n, inSuff)), nx, ny, nz, fmt);
end

zooSum = zeros(nx, ny, nz);
for n = v05zoo
    zooSum = zooSum + read3d(fullfile(inDir, ...
        sprintf('ptracers_optimized_%02d%s', n, inSuff)), nx, ny, nz, fmt);
end

nPhotoV06 = 6;   % nPhoto in code/DARWIN_SIZE.h
nZooV06   = 4;   % nplank - nPhoto

if clipNegatives
    phytoSum = max(phytoSum, 0);
    zooSum   = max(zooSum,   0);
end

phytoEach = phytoSum / nPhotoV06;
zooEach   = zooSum   / nZooV06;

for k = 1:nPhotoV06                       % v06 tracers 21..26 = c01..c06
    n = 20 + k;
    write3d(fullfile(outDir, sprintf('ptracers_v06_3deg_%02d%s', n, outSuff)), ...
            phytoEach, fmt);
    fprintf('%2d  c%02d     phyto  (v05 phyto sum / %d)\n', n, k, nPhotoV06);
end

for k = 1:nZooV06                         % v06 tracers 27..30 = c07..c10
    n = 26 + k;
    write3d(fullfile(outDir, sprintf('ptracers_v06_3deg_%02d%s', n, outSuff)), ...
            zooEach, fmt);
    fprintf('%2d  c%02d     zoo    (v05 zoo sum / %d)\n', n, nPhotoV06+k, nZooV06);
end

%% 31-36: Chl01..Chl06 -- not written on purpose
%
% data.ptracers leaves PTRACERS_initialFile(31:36)=' ' and data.darwin sets
% darwin_chlInitBalanced=T, so Darwin derives Chl from biomass. Writing files
% here would override that and reintroduce the v05 Chl:C ratios.

fprintf('31-36 Chl01..Chl06  not written (darwin_chlInitBalanced=T)\n');

%% conservation check

fprintf('\ntotal plankton C conserved:\n');
fprintf('  v05 sum  = %.6e\n', sum(phytoSum(:)) + sum(zooSum(:)));
fprintf('  v06 sum  = %.6e\n', nPhotoV06*sum(phytoEach(:)) + nZooV06*sum(zooEach(:)));

%% ------------------------------------------------------------------------

function fld = read3d(fname, nx, ny, nz, fmt)
    fid = fopen(fname, 'r', fmt);
    if fid < 0
        error('cannot open %s', fname);
    end
    fld = fread(fid, nx*ny*nz, 'real*4=>double');
    fclose(fid);
    if numel(fld) ~= nx*ny*nz
        error('%s: expected %d values, got %d', fname, nx*ny*nz, numel(fld));
    end
    fld = reshape(fld, [nx ny nz]);
end

function write3d(fname, fld, fmt)
    fid = fopen(fname, 'w', fmt);
    if fid < 0
        error('cannot write %s', fname);
    end
    fwrite(fid, fld, 'real*4');
    fclose(fid);
end

function s = tracerName(n)
    names = {'DIC','NO3','NO2','NH4','PO4','FeT','SiO2','DOC','DON','DOP', ...
             'DOFe','POC','PON','POP','POFe','POSi','PIC','ALK','O2'};
    s = names{n};
end
