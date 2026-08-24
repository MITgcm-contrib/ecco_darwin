%% regrid_loac_runoff_to_3deg.m
%
% TEMPLATE -- regrid the 11 LOAC riverine BGC forcing fields from the ECCO
% V4r5 (llc90) grid onto the 128x64 3deg grid, so that DARWIN_NUTRIENT_RUNOFF
% can be enabled in v06/3deg.
%
% This is a template, not a turnkey script: it needs gcmfaces on the MATLAB
% path and the llc90 grid directory, and the conservative-remap step is the
% part you should check before trusting the output. Runoff is a coastal,
% highly localized field -- a naive bilinear interpolation from llc90 to a
% 2.8125-degree grid will smear river mouths across whole basins and will NOT
% conserve the integrated flux. Use area-weighted binning (below) instead.
%
% After running this, in v06/3deg:
%   1. code/DARWIN_OPTIONS.h  -> #define DARWIN_NUTRIENT_RUNOFF
%   2. input/data.darwin      -> uncomment the runoff block and the
%                                *runoff_interpMethod lines
%   3. rebuild
%
% Source (NAS): /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/

%% settings

srcDir  = '/nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/';
gridDir = '/nobackup/dcarrol2/grid/V4r5/';           % llc90 grid for gcmfaces
outDir  = '/nobackup/dcarrol2/pub/3deg/v06/bgc_runoff/';

% the 11 fields referenced by input/data.darwin
flds = {'DIC_ECCO_V4r5','DIN_ECCO_V4r5','DIP_ECCO_V4r5','DSi_ECCO_V4r5', ...
        'DOC_ECCO_V4r5','DON_ECCO_V4r5','DOP_ECCO_V4r5', ...
        'POC_ECCO_V4r5','PON_ECCO_V4r5','PP_ECCO_V4r5'};
% note: data.darwin lists POPrunofffile = 'PP_ECCO_V4r5' (PP is the POP field)

% target grid: 128x64 global lat-lon, matching code/SIZE.h and input/data
%   delX = 128*2.8125, delY = 64*2.8125, ygOrigin = -90
nx = 128;
ny =  64;
dx = 2.8125;
lonEdge = -180 : dx : 180;          % 129 edges
latEdge =  -90 : dx : 90;           %  65 edges
lonC    = lonEdge(1:end-1) + dx/2;
latC    = latEdge(1:end-1) + dx/2;

nDays = 12053;                      % 1992-01-01 .. 2024-12-31, daily records
fmt   = 'ieee-be';

if ~exist(outDir,'dir'); mkdir(outDir); end

%% load the source grid

global mygrid
if isempty(mygrid)
    grid_load(gridDir, 5, 'compact');
end

srcLon  = convert2gcmfaces(mygrid.XC);
srcLat  = convert2gcmfaces(mygrid.YC);
srcArea = convert2gcmfaces(mygrid.RAC);

srcLon  = srcLon(:);
srcLat  = srcLat(:);
srcArea = srcArea(:);

% precompute which target cell each source cell falls into -- runoff is a
% flux per unit area, so we bin flux*area and divide by target area to
% conserve the global integral exactly.
ix = discretize(srcLon, lonEdge);
iy = discretize(srcLat, latEdge);
good = ~isnan(ix) & ~isnan(iy);
binIdx = sub2ind([nx ny], ix(good), iy(good));

% target cell areas on a sphere
Re = 6371000;
dLon = dx * pi/180;
tgtArea = zeros(nx, ny);
for j = 1:ny
    tgtArea(:,j) = Re^2 * dLon * (sind(latEdge(j+1)) - sind(latEdge(j)));
end

%% regrid each field, one daily record at a time

for f = 1:numel(flds)
    src = fullfile(srcDir, flds{f});
    dst = fullfile(outDir, strrep(flds{f}, 'ECCO_V4r5', '3deg'));

    fidI = fopen(src, 'r', fmt);
    fidO = fopen(dst, 'w', fmt);
    if fidI < 0, error('cannot open %s', src); end

    nSrc = numel(srcArea);
    totIn = 0; totOut = 0;

    for d = 1:nDays
        rec = fread(fidI, nSrc, 'real*4=>double');
        if numel(rec) < nSrc
            fprintf('%s: stopped after %d records\n', flds{f}, d-1);
            break
        end

        flux = rec .* srcArea;                       % per-area -> total
        acc  = accumarray(binIdx, flux(good), [nx*ny 1], @sum, 0);
        out  = reshape(acc, [nx ny]) ./ tgtArea;     % back to per-area

        fwrite(fidO, out, 'real*4');

        totIn  = totIn  + sum(flux(good));
        totOut = totOut + sum(acc);
    end

    fclose(fidI);
    fclose(fidO);

    fprintf('%-16s  flux conserved to %.3e relative\n', ...
            flds{f}, abs(totOut-totIn)/max(abs(totIn),eps));
end

%% notes
%
% * The binning above conserves the global integral but does nothing about
%   land/ocean mismatch: an llc90 coastal cell can land in a 3deg cell that
%   the 3deg bathymetry calls land, and that flux is then lost when the model
%   masks it. Check the loss by comparing sum(out .* oceanMask3deg .* tgtArea)
%   against totIn, and if it is large, nudge stranded flux to the nearest wet
%   3deg cell before writing.
%
% * At 2.8125 degrees most individual rivers are unresolved. Whether runoff
%   should be on at all in a 3deg configuration is a judgment call -- the
%   default v06/3deg setup leaves it off for exactly this reason.
