function dyfamed_mat2csv(site_dir)
% Convert DYFAMED raw files that Python cannot read directly into plain CSV.
%   - SEANOE 89000.mat (MATLAB table object, bottle data 1994-2019)
%   - SEANOE DYFAMED_trap_fluxes_1988_2017.xls (inside 134464.tgz)
%   - obs-vlfr SO-DYFAMED .xls files (pigments, PP, POC/PON, nutrients)
% Usage (from MATLAB, batch):
%   matlab -batch "addpath('<obs>/scripts'); dyfamed_mat2csv('<obs>/DYFAMED')"
% Only data are loaded (load/readcell); nothing in raw/ is executed.
% Output: <site_dir>/raw/<source>/csv_matlab/*.csv (read by std_dyfamed.py)
raw = fullfile(site_dir, 'raw');

% 1. SEANOE bottle table
out = fullfile(raw, 'seanoe', 'csv_matlab');
if ~exist(out, 'dir'), mkdir(out); end
S = load(fullfile(raw, 'seanoe', '89000.mat'), 'DYFAMED_BTLE', 'QCcodes');
writetable(S.DYFAMED_BTLE, fullfile(out, '89000_DYFAMED_BTLE.csv'));
writecell(cellfun(@string, S.QCcodes, 'UniformOutput', false), fullfile(out, '89000_QCcodes.csv'));

% 2. SEANOE trap xls (sheet 1)
xls_all(fullfile(raw, 'seanoe', 'x134464', 'dyfamed-mooring', 'DYFAMED_trap_fluxes_1988_2017.xls'), out);

% 3. obs-vlfr xls files
vin = fullfile(raw, 'obsvlfr_sodyf');
vout = fullfile(vin, 'csv_matlab');
if ~exist(vout, 'dir'), mkdir(vout); end
f = dir(fullfile(vin, '*.xls'));
for k = 1:numel(f)
    xls_all(fullfile(vin, f(k).name), vout);
end
end

function xls_all(inf, outdir)
[~, base] = fileparts(inf);
sh = sheetnames(inf);
for i = 1:numel(sh)
    C = readcell(inf, 'Sheet', sh{i});
    if isempty(C), continue; end
    m = cellfun(@(x) isa(x, 'missing'), C);
    C(m) = {''};
    d = cellfun(@(x) isdatetime(x), C);
    C(d) = cellfun(@(x) char(string(x, 'yyyy-MM-dd HH:mm:ss')), C(d), 'UniformOutput', false);
    nm = regexprep(sh{i}, '[^A-Za-z0-9_-]', '_');
    writecell(C, fullfile(outdir, sprintf('%s__%02d_%s.csv', base, i, nm)));
end
end
