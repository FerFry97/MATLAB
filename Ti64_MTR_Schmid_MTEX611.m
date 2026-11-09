%% Ti-6Al-4V: MTR (macrozone) detection + Schmid factor of the MTRs, MTEX 6.1.1
% Stand-alone program: run it separately for each map (e.g. 2.5 um and 0.15 um).
%  1. Detects MTRs / macrozones in the .ctf with the method of
%     B. Fernandez Silva et al., J. Appl. Cryst. 56 (2023) 737-749
%     (same code as the CTF MACROZONES tool; MTR map checked against MATLAB results)
%  2. Computes basal / prismatic / pyramidal Schmid factors, c-axis-to-load angle
%     and basal normal stress ratio for every indexed alpha pixel (MTEX)
%  3. Looks at the MTRs in isolation: MTR vs non-MTR, per macrozone, and how each
%     macrozone responds when the load direction is rotated in the map plane.
% Geometric, uniaxial analysis: no local stress field, no slip prediction.
% NOT TESTED in MATLAB (written without MATLAB/MTEX access).
%
% OUTPUT (folder next to the .ctf: <name>_MTR_Schmid)
%  M01 macrozone map          M03-M06 maps of c-axis angle / SF inside macrozones
%  M07 MTR vs non-MTR distributions      M08 load-direction sweep per macrozone
%  T1 summary (all / MTR / macrozones / non-MTR)   T2 per macrozone
%  T3 load sweep (angle x macrozone)                T4 load sensitivity per macrozone
%  MTR_Macrozones_data/statistics.csv, MTR_Schmid_Results.mat

clear; close all; clc;

%% 0. Start MTEX
% This MUST be the full MTEX installation root (contains startup_mtex.m
% AND tools/option_tools/setMTEXpref.m), NOT the folder with your own scripts.
% EDIT mtexRoot if your installation is located somewhere else.
mtexRoot = 'C:\mtex-6.1.1';
if ~isappdata(0,'mtex') || exist('EBSD','class') ~= 8
    startupFile = fullfile(mtexRoot,'startup_mtex.m');
    prefFile = fullfile(mtexRoot,'tools','option_tools','setMTEXpref.m');
    if exist(startupFile,'file') ~= 2 || exist(prefFile,'file') ~= 2
        error(['MTEX installation path is incorrect or incomplete. ' ...
               'Edit mtexRoot at the top of this script so it points to ' ...
               'the MTEX ROOT directory. Current value: %s'],mtexRoot);
    end
    % startup_mtex must run FROM its own folder: it discovers its toolbox
    % subfolders relative to the location of startup_mtex.m.
    originalFolder = pwd;
    try
        cd(mtexRoot);
        startup_mtex('noMenu');
    catch ME
        cd(originalFolder);
        rethrow(ME);
    end
    cd(originalFolder);
end
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
setMTEXpref('defaultColorMap',WhiteJetColorMap);

%% 1. Options
alphaPhase = 'Ti-Hex';               % MTEX phase name of alpha in the .ctf

% --- Load (EBSD map frame; X east, Y down, Z into the plane) ---
opt = struct();
opt.loadAxisXYZ = [-1 0 0];          % uniaxial load direction
opt.cAxisClassDeg = [30 60];         % c-axis-to-load: <30 'hard', 30-60 'intermediate', >60 'soft'
opt.loadSweep = true;                % M08: rotate the load in the map plane
opt.loadSweepStepDeg = 5;            % 0:step:180 deg from +X

% --- MTR detection (Fernandez Silva et al. 2023; paper defaults) ---
opt.MTR.windowUm      = 50;          % window size in um (gridsize = 'auto')
opt.MTR.gridsize      = 'auto';      % or an ODD number of pixels of the detection grid
opt.MTR.critDisor     = 20;          % deg
opt.MTR.critFraction  = 0.3;
opt.MTR.areaThreshold = 10000;       % um^2 (lower it for small fine-step maps)
% Detection grid: 'auto' = native step if >= 1 um; for finer maps (e.g. 0.15 um)
% the MTRs are detected on a resampled grid with this step (every k-th pixel)
% and the result is assigned to every pixel of the fine map.
opt.MTR.detectStepUm  = 'auto';      % 'auto' or a number in um (e.g. 0.5, 1, 2.5)
opt.MTR.autoFineStepUm = 1.0;        % used by 'auto' when the native step < 1 um

%% 2. Select the .ctf file(s)
[f,p] = uigetfile({'*.ctf','Oxford / AZtec CTF (*.ctf)'}, ...
    'Select one or more .ctf maps','MultiSelect','on');
if isequal(f,0), disp('No file selected.'); return; end
files = fullfile(p,cellstr(f));

%% 3. Analysis, one map at a time
for n = 1:numel(files)
    fname = files{n};
    [~,sampleName] = fileparts(fname);
    outputFolder = fullfile(fileparts(fname),[sampleName '_MTR_Schmid']);
    if ~isfolder(outputFolder), mkdir(outputFolder); end
    fprintf('\n================ %s ================\nOUTPUT: %s\n',sampleName,outputFolder);
    try
        Res = analyseMTRSchmid(fname,alphaPhase,opt,outputFolder); %#ok<NASGU>
    catch ME
        warning('Analysis of %s failed: %s',sampleName,ME.message);
        fprintf('%s\n',getReport(ME,'extended','hyperlinks','off'));
    end
end

%% ========================================================================
function Res = analyseMTRSchmid(fname,alphaPhase,opt,folder)
% ---- MTEX data ----
ebsd = EBSD.load(fname,'interface','ctf','convertEuler2SpatialReferenceFrame');
try
    alphaData = ebsd(alphaPhase);
catch
    disp(ebsd); error('Phase "%s" not found: edit alphaPhase.',alphaPhase);
end
if isempty(alphaData), error('No indexed alpha pixels.'); end
CSalpha = alphaData.CS;
fprintf('MTEX: %d indexed alpha pixels\n',length(alphaData));

% ---- MTR detection on the same .ctf ----
G0 = mtr_ReadCTF(fname,'auto');
[G,k] = mtr_DetectionGrid(G0,opt.MTR);
Pm = opt.MTR;
if ischar(Pm.gridsize)
    g = round(Pm.windowUm/G.step); if mod(g,2)==0, g = g+1; end
    Pm.gridsize = max(3,g);
end
mapAreaUm2 = (G0.nx*G0.step)*(G0.ny*G0.stepY);
fprintf('MTR detection grid: step %.3g um (every %d px) | window %d px (%.1f um), %g deg, %g %%, area >= %g um2\n', ...
    G.step,k,Pm.gridsize,Pm.gridsize*G.step,Pm.critDisor,100*Pm.critFraction,Pm.areaThreshold);
if Pm.areaThreshold > 0.25*mapAreaUm2
    warning(['Area threshold (%g um2) is more than 25%% of the map area (%.0f um2): few or no ' ...
        'macrozones can pass it. For small fine-step maps lower opt.MTR.areaThreshold.'], ...
        Pm.areaThreshold,mapAreaUm2);
end
if Pm.gridsize*G.step > 0.5*min(G0.nx*G0.step,G0.ny*G0.stepY)
    warning('The %.0f um window is more than half of the shorter map side.',Pm.gridsize*G.step);
end
M = mtr_DetectMacrozones(G,Pm);
writetable(M.DataTable,fullfile(folder,'MTR_Macrozones_data.csv'));
writetable(M.Statistics,fullfile(folder,'MTR_Macrozones_statistics.csv'),'WriteRowNames',true);
try
    figure('Name','M01 macrozones','Color','white');
    imshow(M.MacroNumberColor_count,[]);
    title(sprintf('M01 macrozones >= %g um^2: %d | MTR %.1f %% of map (detection step %.3g um)', ...
        Pm.areaThreshold,numel(M.kept),M.MTR,G.step),'Interpreter','none');
    saveFigure(folder,'M01_Macrozones_map.png');
catch
end

% ---- every alpha pixel -> MTR / macrozone ----
MZ = mtr_MapAlphaPixels(alphaData,G0,G,k,M,fname);
MZ.areaThreshold = Pm.areaThreshold;

% ---- Schmid factors (same definitions as the reconstruction script) ----
xyz = opt.loadAxisXYZ(:).'; xyz = xyz/norm(xyz);
loadDir = vector3d(xyz(1),xyz(2),xyz(3));
ori = alphaData.orientations;
loadInCrystal = inv(ori) * loadDir;
sBas = slipSystem.basal(CSalpha).symmetrise('antipodal');
sPri = slipSystem.prismaticA(CSalpha).symmetrise('antipodal');
sPyr = slipSystem.pyramidalCA(CSalpha).symmetrise('antipodal');
mBas = max(abs(sBas.SchmidFactor(loadInCrystal)),[],2); mBas = mBas(:);
mPri = max(abs(sPri.SchmidFactor(loadInCrystal)),[],2); mPri = mPri(:);
mPyr = max(abs(sPyr.SchmidFactor(loadInCrystal)),[],2); mPyr = mPyr(:);
cAxisSpec = ori * Miller(0,0,0,1,CSalpha,'UVTW');
theta = angle(cAxisSpec,loadDir)./degree; theta = theta(:);
theta = min(theta,180-theta);
basalNormalRatio = cosd(theta).^2;
mEasyA = max(mBas,mPri);

% ---- T1 summary ----
labels = {'ALL measured alpha','MTR pixels', ...
    sprintf('Macrozones >= %g um2',MZ.areaThreshold),'Alpha outside MTR'};
masks = {true(numel(mBas),1), MZ.inMTR, MZ.id>0, ~MZ.inMTR};
rows = cell(numel(labels),1);
for j = 1:numel(labels)
    rows{j} = schmidSummaryRow(labels{j},masks{j},mBas,mPri,mPyr,mEasyA,theta);
end
T1 = vertcat(rows{:});
writetable(T1,fullfile(folder,'T1_Schmid_Summary_MTR.csv'));
fprintf('\n--- T1 SCHMID FACTOR: ALL / MTR / MACROZONES / NON-MTR (load [%g %g %g]) ---\n',xyz);
disp(T1);

% ---- T2 per macrozone + maps inside macrozones ----
T2 = table();
inMZ = MZ.id > 0;
if any(inMZ)
    T2 = schmidPerMacrozone(MZ,mBas,mPri,mPyr,mEasyA,theta,basalNormalRatio,opt);
    writetable(T2,fullfile(folder,'T2_Schmid_per_Macrozone.csv'));
    fprintf('\n--- T2 SCHMID FACTOR PER MACROZONE ---\n');
    disp(T2(:,{'Macrozone','Area_um2','Basal_mean','Prismatic_mean', ...
        'MaxA_mean','CAxisToLoad_mean_deg','CAxisClass'}));
    a = sprintf(' >= %g um2',MZ.areaThreshold);
    plotAlphaSchmidMap(alphaData(inMZ),theta(inMZ),[0 90],'angle [deg]', ...
        ['M03 macrozones' a ': c-axis to load angle'],'M03_CAxis_Angle_Macrozones.png',xyz,folder);
    plotAlphaSchmidMap(alphaData(inMZ),mEasyA(inMZ),[0 0.5],'m (absolute)', ...
        ['M04 macrozones' a ': max(basal, prismatic) SF'],'M04_SF_MaxA_Macrozones.png',xyz,folder);
    plotAlphaSchmidMap(alphaData(inMZ),mBas(inMZ),[0 0.5],'m (absolute)', ...
        ['M05 macrozones' a ': basal <a> SF'],'M05_SF_Basal_Macrozones.png',xyz,folder);
    plotAlphaSchmidMap(alphaData(inMZ),mPri(inMZ),[0 0.5],'m (absolute)', ...
        ['M06 macrozones' a ': prismatic <a> SF'],'M06_SF_Prismatic_Macrozones.png',xyz,folder);
else
    warning('No macrozone passed the area threshold: only T1 and M07 are produced.');
end

% ---- M07 MTR vs non-MTR, M08 load-direction sweep ----
try
    plotMTRvsNonMTR(MZ,mBas,mPri,theta,folder);
catch MEh
    warning('M07 failed: %s',MEh.message);
end
S = [];
if opt.loadSweep && any(inMZ)
    try
        S = loadSweepMacrozones(ori,sBas,sPri,MZ,xyz,opt,folder,CSalpha);
    catch MEs
        warning('M08 failed: %s',MEs.message);
    end
end

% ---- save ----
Res = struct('file',fname,'opt',opt,'MTRparams',Pm,'detectStep_um',G.step, ...
    'summary',T1,'perMacrozone',T2,'loadSweep',S,'mapping',rmfield(MZ,'dataTable'), ...
    'x_um',alphaData.x(:),'y_um',alphaData.y(:),'SF_Basal',mBas,'SF_Prismatic',mPri, ...
    'SF_PyramidalCA',mPyr,'CAxisToLoad_deg',theta,'BasalNormalStressOverSigma',basalNormalRatio);
M.MacroNumberColor = []; %#ok<STRNU>
save(fullfile(folder,'MTR_Schmid_Results.mat'),'Res','M','-v7.3');
fprintf('Saved in: %s\n',folder);
end

%% ---------- detection grid (resampling for fine maps) ---------------------
function [G,k] = mtr_DetectionGrid(G0,P)
% k = 1: native grid. k > 1: every k-th pixel (fine maps), step = k*native step.
if ischar(P.detectStepUm)
    if G0.step < 1, target = P.autoFineStepUm; else, target = G0.step; end
else
    target = P.detectStepUm;
end
k = max(1,round(target/G0.step));
G = G0;
if k == 1, return; end
r = 1:k:G0.ny; c = 1:k:G0.nx;
[C,Rr] = meshgrid(c,r);
lin = sub2ind([G0.ny G0.nx],Rr(:),C(:));
G.ny = numel(r); G.nx = numel(c);
G.step = G0.step*k; G.stepY = G0.stepY*k;
G.phase = G0.phase(r,c);
G.euler = G0.euler(r,c,:);
G.indexed = G0.indexed(r,c);
G.G = G0.G(lin,:);
G.Q = G0.Q(lin,:);
fprintf('Fine map: MTR detection on every %d-th pixel (%d x %d grid, step %.3g um)\n',k,G.nx,G.ny,G.step);
end

%% ---------- alpha pixel -> MTR / macrozone --------------------------------
function MZ = mtr_MapAlphaPixels(alphaData,G0,G,k,M,fname)
% The native grid (G0) of the same .ctf fixes the mapping: the orientation
% (coordinate sign / axis flips) is the one where alpha pixels land on Ti-Hex.
[minX,minY] = ctfMinXY(fname);
x = alphaData.x(:); y = alphaData.y(:);
cand = {}; names = {};
for sx = [1 -1]
    for sy = [1 -1]
        c0 = round((sx*x - minX)/G0.step)  + 1;
        r0 = round((sy*y - minY)/G0.stepY) + 1;
        base = sprintf('x%s y%s',char('+'*(sx>0)+'-'*(sx<0)),char('+'*(sy>0)+'-'*(sy<0)));
        cand(end+1,:) = {r0, c0};                  names{end+1} = base; %#ok<AGROW>
        cand(end+1,:) = {G0.ny+1-r0, c0};          names{end+1} = [base ', flip rows']; %#ok<AGROW>
        cand(end+1,:) = {r0, G0.nx+1-c0};          names{end+1} = [base ', flip cols']; %#ok<AGROW>
        cand(end+1,:) = {G0.ny+1-r0, G0.nx+1-c0};  names{end+1} = [base ', flip both']; %#ok<AGROW>
    end
end
score = zeros(size(cand,1),1);
for i = 1:size(cand,1)
    rr = cand{i,1}; cc = cand{i,2};
    ok = rr>=1 & rr<=G0.ny & cc>=1 & cc<=G0.nx;
    if ~any(ok), continue; end
    score(i) = mean(G0.indexed(sub2ind([G0.ny G0.nx],rr(ok),cc(ok)))) * mean(ok);
end
[best,i] = max(score);
fprintf('MTEX <-> MTR grid: %.1f%% of alpha pixels on Ti-Hex (mapping: %s)\n',100*best,names{i});
if best < 0.95
    error('MTEX and MTR grids do not match (%.1f%%). Same .ctf?',100*best);
end
rr = cand{i,1}; cc = cand{i,2};
ok = rr>=1 & rr<=G0.ny & cc>=1 & cc<=G0.nx;
% native pixel -> nearest pixel of the (possibly resampled) detection grid
rd = min(G.ny, round((rr(ok)-1)/k) + 1);
cd = min(G.nx, round((cc(ok)-1)/k) + 1);
lin = sub2ind([G.ny G.nx],rd,cd);
n = numel(x);
MZ.inMTR = false(n,1); MZ.inMTR(ok) = M.ismacro(lin);
MZ.id = zeros(n,1);    MZ.id(ok) = M.MacroNumber_count(lin);
MZ.dataTable = M.DataTable;
MZ.gridMatch = best; MZ.flip = names{i};
fprintf('  alpha pixels in MTR: %d (%.1f%%) | in macrozones: %d (%.1f%%)\n', ...
    nnz(MZ.inMTR),100*mean(MZ.inMTR),nnz(MZ.id>0),100*mean(MZ.id>0));
end

%% ---------- plotting / tables (same as the reconstruction script v5.2) -----
function plotAlphaSchmidMap(e,mapVals,lim,colorbarTitle,ttl,filename,xyz,folder)
try
    figure('Name',ttl,'Color','white');
    plot(e,mapVals);
    % Dark background for areas without indexed alpha
    set(gca,'Color',[0.12 0.12 0.12]); 
    colormap(gca,parula(256)); setColorRange(lim);
    mtexColorbar('title',colorbarTitle);
    tt = title(sprintf('%s | load axis [%g %g %g]',ttl,xyz),'Interpreter','none');
    % Load arrow in a thin band BETWEEN the title and the map (never on the map).
    % Illustrative only, not a stress-vector simulation.
    drawLoadArrowAboveMap(gca,tt,xyz);
    saveFigure(folder,filename);
catch ME
    warning('Could not save %s: %s',filename,ME.message);
end
end

function drawLoadArrowAboveMap(ax,tt,xyz)
% Opens a small gap between the map and its title and draws the load arrow
% in it, centred over the map. Screen directions follow the MTEX settings of
% this script (xAxisDirection east, zAxisDirection intoPlane => +Y points DOWN).
d = xyz(1:2);
if norm(d) < 1e-10, return; end          % load along Z: no in-plane arrow
d = d/norm(d);
gap = 0.06;                               % band height (normalized figure units)
try
    drawnow;
    fig = ancestor(ax,'figure');
    set(ax,'Units','normalized');
    p = get(ax,'Position');
    % make room: lower the top of the axes by the band height
    set(ax,'Position',[p(1) p(2) p(3) max(0.2,p(4)-gap)]);
    drawnow;
    pb = plotBoxNormalized(ax);           % visible map box [x y w h]
    p  = get(ax,'Position');
    % title just above the band
    set(tt,'Units','normalized');
    tp = get(tt,'Position');
    tp(2) = (pb(2) + pb(4) + gap - p(2)) / p(4);
    set(tt,'Position',tp,'VerticalAlignment','bottom');
    % arrow centred over the map, in the middle of the band
    yc = pb(2) + pb(4) + 0.5*gap;
    xc = pb(1) + 0.5*pb(3);
    L  = 0.10;                            % arrow length (fraction of figure width)
    dx = L*d(1);  dy = -0.5*gap*d(2);     % screen up = -Y
    xa = [xc - dx/2, xc + dx/2];  ya = [yc - dy/2, yc + dy/2];
    annotation(fig,'arrow',xa,ya,'Color','black','LineWidth',1.6, ...
        'HeadLength',8,'HeadWidth',8);
    annotation(fig,'textbox',[max(xa)+0.008, yc-0.02, 0.08, 0.04], ...
        'String','Load','EdgeColor','none','FontSize',10, ...
        'VerticalAlignment','middle','HorizontalAlignment','left', ...
        'FitBoxToText','off');
catch MEarrow
    warning('Load-arrow annotation failed: %s',MEarrow.message);
end
end

function pb = plotBoxNormalized(ax)
% Position of the visible plot box (after axis equal/image) in normalized
% figure units. Falls back to the axes Position if the aspect is free.
set(ax,'Units','pixels'); pos = get(ax,'Position'); set(ax,'Units','normalized');
fig = ancestor(ax,'figure');
oldU = get(fig,'Units'); set(fig,'Units','pixels'); fp = get(fig,'Position'); set(fig,'Units',oldU);
pb = pos;
if strcmp(get(ax,'DataAspectRatioMode'),'manual') || strcmp(get(ax,'PlotBoxAspectRatioMode'),'manual')
    xl = diff(get(ax,'XLim')); yl = diff(get(ax,'YLim'));
    dar = get(ax,'DataAspectRatio');
    w = xl/dar(1); h = yl/dar(2);
    if strcmp(get(ax,'PlotBoxAspectRatioMode'),'manual') && ~strcmp(get(ax,'DataAspectRatioMode'),'manual')
        pbar = get(ax,'PlotBoxAspectRatio'); w = pbar(1); h = pbar(2);
    end
    ar = w/h;
    if pos(3)/pos(4) > ar           % limited by height
        nw = pos(4)*ar; pb = [pos(1)+(pos(3)-nw)/2, pos(2), nw, pos(4)];
    else                            % limited by width
        nh = pos(3)/ar; pb = [pos(1), pos(2)+(pos(4)-nh)/2, pos(3), nh];
    end
end
pb = pb ./ [fp(3) fp(4) fp(3) fp(4)];
end

function saveFigure(folder,name)
exportgraphics(gcf,fullfile(folder,name),'Resolution',300);
end

function [minX,minY] = ctfMinXY(fname)
% Minimum X and Y of all points of a .ctf (grid origin used by the macrozone tool).
fid = fopen(fname,'r');
if fid < 0, error('Cannot open %s',fname); end
c = onCleanup(@() fclose(fid));
cols = {};
while true
    ln = fgetl(fid);
    if ~ischar(ln), error('Column line "Phase X Y ..." not found in %s',fname); end
    if ~isempty(regexp(ln,'^Phase\t','once'))
        cols = strsplit(strtrim(ln),'\t'); break
    end
end
A = fscanf(fid,'%f',[numel(cols) Inf]);
minX = min(A(strcmpi(cols,'X'),:));
minY = min(A(strcmpi(cols,'Y'),:));
end

function row = schmidSummaryRow(label,mask,b,p,c,a,theta)
if ~any(mask)
    row = table({label},0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
        'VariableNames',{'Subset','Pixels','Basal_mean','Basal_P90', ...
        'Prismatic_mean','Prismatic_P90','PyramidalCA_mean', ...
        'PyramidalCA_P90','MaxA_mean','CAxisAngle_mean_deg'});
else
    row = table({label},nnz(mask),mean(b(mask)),quantile(b(mask),0.9), ...
        mean(p(mask)),quantile(p(mask),0.9),mean(c(mask)),quantile(c(mask),0.9), ...
        mean(a(mask)),mean(theta(mask)), ...
        'VariableNames',{'Subset','Pixels','Basal_mean','Basal_P90', ...
        'Prismatic_mean','Prismatic_P90','PyramidalCA_mean', ...
        'PyramidalCA_P90','MaxA_mean','CAxisAngle_mean_deg'});
end
end

function T = schmidPerMacrozone(MZ,b,p,c,a,theta,sn,opt)
% One row per macrozone: mean and P90 Schmid factors, c-axis angle, class.
ids = unique(MZ.id(MZ.id>0));
n = numel(ids);
z = nan(n,1);
T = table(ids,z,z,z,z,z,z,z,z,z,z,z,z,repmat({''},n,1), ...
    'VariableNames',{'Macrozone','AlphaPixels','Area_um2','Basal_mean','Basal_P90', ...
    'Prismatic_mean','Prismatic_P90','PyramidalCA_mean','MaxA_mean','MaxA_P90', ...
    'CAxisToLoad_mean_deg','CAxisToLoad_std_deg','BasalSigmaN_over_sigma_mean','CAxisClass'});
cl = opt.cAxisClassDeg;
DT = MZ.dataTable;
for i = 1:n
    m = MZ.id == ids(i);
    T.AlphaPixels(i) = nnz(m);
    j = find(DT.Macrozone == ids(i),1);
    if ~isempty(j), T.Area_um2(i) = DT.Area_um2(j); end
    T.Basal_mean(i) = mean(b(m));      T.Basal_P90(i) = quantile(b(m),0.9);
    T.Prismatic_mean(i) = mean(p(m));  T.Prismatic_P90(i) = quantile(p(m),0.9);
    T.PyramidalCA_mean(i) = mean(c(m));
    T.MaxA_mean(i) = mean(a(m));       T.MaxA_P90(i) = quantile(a(m),0.9);
    T.CAxisToLoad_mean_deg(i) = mean(theta(m));
    T.CAxisToLoad_std_deg(i) = std(theta(m));
    T.BasalSigmaN_over_sigma_mean(i) = mean(sn(m));
    if T.CAxisToLoad_mean_deg(i) < cl(1)
        T.CAxisClass{i} = 'hard (c ~ parallel to load)';
    elseif T.CAxisToLoad_mean_deg(i) > cl(2)
        T.CAxisClass{i} = 'soft (c ~ perpendicular to load)';
    else
        T.CAxisClass{i} = 'intermediate';
    end
end
T = sortrows(T,'Area_um2','descend','MissingPlacement','last');
end

function plotMTRvsNonMTR(MZ,b,p,theta,folder)
% F16: distributions of basal / prismatic SF and c-axis angle, MTR vs non-MTR.
inM = MZ.inMTR; inZ = MZ.id > 0;
groups = {~inM, inM, inZ};
names = {'alpha outside MTR','MTR pixels',sprintf('macrozones >= %g um2',MZ.areaThreshold)};
vals = {b, p, theta}; ttl = {'Basal <a> m','Prismatic <a> m','c-axis to load [deg]'};
edges = {linspace(0,0.5,26), linspace(0,0.5,26), linspace(0,90,19)};
figure('Name','M07 MTR vs non-MTR','Color','white','Position',[80 80 1350 420]);
for k = 1:3
    subplot(1,3,k); hold on
    for g = 1:3
        if ~any(groups{g}), continue; end
        h = histcounts(vals{k}(groups{g}),edges{k},'Normalization','probability');
        c = 0.5*(edges{k}(1:end-1)+edges{k}(2:end));
        plot(c,h,'LineWidth',1.8);
    end
    hold off; grid on; xlabel(ttl{k}); ylabel('fraction of pixels');
    if k == 1, legend(names,'Location','best'); end
    title(ttl{k});
end
sgtitle('M07 Schmid factor of MTRs in isolation vs the rest of alpha (same load axis)');
saveFigure(folder,'M07_SF_MTR_vs_nonMTR.png');
end

function S = loadSweepMacrozones(ori,sBas,sPri,MZ,xyz,opt,folder,CSalpha)
% M08: mean basal / prismatic SF and c-axis angle of every macrozone as the
% uniaxial load direction is rotated in the map plane (0-180 deg from +X).
% Shows how sensitive each MTR is to the loading direction.
ids = unique(MZ.id(MZ.id>0));
ang = 0:opt.loadSweepStepDeg:180;
nA = numel(ang); nZ = numel(ids);
Bm = zeros(nA,nZ); Pm = Bm; Am = Bm; Th = Bm;
cDir = Miller(0,0,0,1,CSalpha,'UVTW');
for z = 1:nZ
    m = MZ.id == ids(z);
    o = ori(m);
    cAx = o * cDir;
    for i = 1:nA
        d = vector3d(cosd(ang(i)),sind(ang(i)),0);
        lc = inv(o) * d;
        mb = max(abs(sBas.SchmidFactor(lc)),[],2);
        mp = max(abs(sPri.SchmidFactor(lc)),[],2);
        t = angle(cAx,d)./degree; t = min(t,180-t);
        Bm(i,z) = mean(mb); Pm(i,z) = mean(mp); Am(i,z) = mean(max(mb,mp)); Th(i,z) = mean(t);
    end
end
a0 = mod(atan2d(xyz(2),xyz(1)),180);   % actual load direction (|m| is antipodal)
lab = arrayfun(@(k) sprintf('MZ %d',k),ids,'UniformOutput',false);
figure('Name','M08 load-direction sweep','Color','white','Position',[80 80 1350 420]);
data = {Bm, Pm, Th}; yl = {'mean basal <a> m','mean prismatic <a> m','mean c-axis to load [deg]'};
lim = {[0 0.5],[0 0.5],[0 90]};
for k = 1:3
    subplot(1,3,k);
    plot(ang,data{k},'LineWidth',1.5); hold on
    xline(a0,'k--','actual load','LabelVerticalAlignment','bottom');
    hold off; grid on; xlim([0 180]); ylim(lim{k});
    xlabel('in-plane load direction [deg from +X]'); ylabel(yl{k});
    if k == 1, legend(lab,'Location','best','FontSize',8); end
end
sgtitle('M08 Sensitivity of each macrozone to the in-plane load direction (uniaxial, geometric)');
saveFigure(folder,'M08_LoadSweep_Macrozones.png');
% table: one row per angle, columns per macrozone
T = table(ang(:),'VariableNames',{'LoadAngle_deg'});
for z = 1:nZ
    T.(sprintf('MZ%d_Basal',ids(z))) = Bm(:,z);
    T.(sprintf('MZ%d_Prism',ids(z))) = Pm(:,z);
    T.(sprintf('MZ%d_MaxA',ids(z)))  = Am(:,z);
    T.(sprintf('MZ%d_CAxis_deg',ids(z))) = Th(:,z);
end
writetable(T,fullfile(folder,'T3_LoadSweep_Macrozones.csv'));
% per macrozone: softest / hardest in-plane direction
[mx,ix] = max(Am,[],1); [mn,in] = min(Am,[],1);
[~,i0] = min(abs(ang - a0));
R = table(ids,Am(i0,:).',mx.',ang(ix).',mn.',ang(in).',max(Am,[],1).'-min(Am,[],1).', ...
    'VariableNames',{'Macrozone','MaxA_at_actual_load','MaxA_max','Angle_of_max_deg', ...
    'MaxA_min','Angle_of_min_deg','MaxA_range'});
writetable(R,fullfile(folder,'T4_LoadSensitivity_Macrozones.csv'));
fprintf('\n--- LOAD-DIRECTION SENSITIVITY PER MACROZONE (M08, T3/T4) ---\n');
disp(R);
S = struct('angle_deg',ang,'ids',ids,'basal',Bm,'prismatic',Pm,'maxA',Am,'cAxis_deg',Th,'summary',R);
end


%% ===================== BUILT-IN MTR DETECTION ============================
% Copied from the CTF MACROZONES tool (RUN_MACROZONES), renamed mtr_* so it
% never clashes with other functions on the MATLAB path. Same algorithm,
% checked against the published results of Fernandez Silva et al. (2023).
function ebsd = mtr_ReadCTF(filename, analysisPhases)
%READ_CTF_MACROZONES  Reads an Oxford .ctf file (AZtec / Channel 5 or
%   exported from MTEX) directly, without Import_CTF or an intermediate .mat.
%
%   ebsd = mtr_ReadCTF(filename)
%   ebsd = mtr_ReadCTF(filename, analysisPhases)
%
%   filename        full path to the .ctf file
%   analysisPhases  'auto' (default) -> use the hexagonal phase(s) (gamma = 120 deg)
%                   or a vector with phase numbers from the .ctf, e.g. 1 or [1 2]
%
%   Differences with the original Import_CTF.m:
%   - The header is parsed by keyword (XCells, XStep, Phases...) instead of
%     by line number, so headers written by MTEX also work.
%   - Every point is placed on the grid from its X,Y coordinates, so cropped
%     maps or maps with missing points work.
%   - Only pixels of the hexagonal phase(s) (alpha-Ti) are treated as
%     "indexed" for the macrozone analysis. Beta (cubic) pixels are treated as
%     non-indexed (the original code analysed them as if they were hexagonal).
%
%   Output (struct):
%     ebsd.name, ebsd.ny, ebsd.nx, ebsd.step, ebsd.stepY
%     ebsd.phase     ny x nx  (0 = non-indexed / no data)
%     ebsd.euler     ny x nx x 3 (degrees, Bunge)
%     ebsd.indexed   ny x nx logical (phase used in the analysis)
%     ebsd.G         (ny*nx) x 9 flattened orientation matrices (as Gmatrix.m)
%     ebsd.Q         (ny*nx) x 4 quaternions (as Quaternion_Euler_Deg.m)
%     ebsd.phases    struct with name, lattice parameters and isHex flag
%     ebsd.hexPhase  phase numbers used in the analysis

if nargin < 1 || isempty(filename)
    [f, d] = uigetfile('*.ctf', 'Select a .ctf file');
    if isequal(f, 0), ebsd = []; disp('No file selected.'); return, end
    filename = fullfile(d, f);
end
if nargin < 2 || isempty(analysisPhases)
    analysisPhases = 'auto';
end
if ~exist(filename, 'file')
    error('File not found: %s', filename);
end

[~, name, ~] = fileparts(filename);
txt = fileread(filename);
txt = strrep(txt, char(13), '');           % remove CR (Windows line endings)

%% Header: find the column title line "Phase X Y ..."
[i0, i1] = regexp(txt, '(^|\n)Phase\t[^\n]*', 'once');
if isempty(i0)
    error('Column line "Phase  X  Y ..." not found in %s', filename);
end
header  = txt(1:i0);
colLine = strtrim(txt(i0:i1));
body    = txt(i1+1:end);
lines = regexp(header, '\n', 'split');
lines = lines(~cellfun(@isempty, strtrim(lines)));

XStep = mtr_readValue(lines, 'XStep');
YStep = mtr_readValue(lines, 'YStep');
jobMode = '';
for k = 1:numel(lines)
    t = regexp(lines{k}, '^JobMode\s+(\S+)', 'tokens', 'once');
    if ~isempty(t), jobMode = t{1}; end
end
if ~isempty(jobMode) && ~strcmpi(jobMode, 'Grid')
    error(['The map uses JobMode = %s. This tool needs a square grid ' ...
           '(JobMode Grid). Convert the hexagonal grid first ' ...
           '(e.g. in MTEX with ebsd.gridify).'], jobMode);
end

%% Phases
iPh = find(strncmp(lines, 'Phases', 6), 1);
nPh = 0;
if ~isempty(iPh)
    nPh = sscanf(lines{iPh}(7:end), '%d');
end
phases = struct('id', {}, 'name', {}, 'abc', {}, 'ang', {}, 'isHex', {});
for k = 1:nPh
    parts = regexp(lines{iPh + k}, '\t', 'split');
    abc = sscanf(strrep(parts{1}, ';', ' '), '%f')';
    ang = sscanf(strrep(parts{2}, ';', ' '), '%f')';
    nm = '';
    if numel(parts) >= 3, nm = strtrim(parts{3}); end
    isHex = numel(ang) == 3 && abs(ang(3) - 120) < 1 && abs(abc(1) - abc(2)) < 1e-3;
    phases(k).id = k; phases(k).name = nm; phases(k).abc = abc;
    phases(k).ang = ang; phases(k).isHex = isHex;
end

if ischar(analysisPhases) && strcmpi(analysisPhases, 'auto')
    hexPhase = [phases([phases.isHex]).id];
    if isempty(hexPhase)
        error(['No hexagonal phase found in the header. ' ...
               'Set the phase number manually in analysisPhases.']);
    end
else
    hexPhase = analysisPhases(:)';
end

%% Data
cols = regexp(colLine, '\t', 'split');
nCol = numel(cols);
data = sscanf(body, '%f');
nRow = floor(numel(data) / nCol);
if nRow * nCol ~= numel(data)
    warning('Number of values is not a multiple of %d columns; extra values ignored.', nCol);
end
data = reshape(data(1:nRow*nCol), nCol, nRow)';

col = @(n) find(strcmpi(cols, n), 1);
cP = col('Phase'); cX = col('X'); cY = col('Y');
cE = [col('Euler1'), col('Euler2'), col('Euler3')];
if isempty(cP) || isempty(cX) || isempty(cY) || numel(cE) < 3
    error('Missing Phase/X/Y/Euler1-3 columns in the .ctf');
end

X = data(:, cX); Y = data(:, cY);
ph = data(:, cP);
E = data(:, cE);

% Grid steps (estimated from coordinates if missing in the header)
if isnan(XStep) || XStep <= 0, XStep = mtr_minStep(X); end
if isnan(YStep) || YStep <= 0, YStep = mtr_minStep(Y); end
if abs(XStep - YStep) > 1e-3 * XStep
    warning('XStep (%.4g) and YStep (%.4g) differ; XStep is used for areas.', XStep, YStep);
end

% Euler angles in radians (some exports) -> degrees
iE = ph > 0;
if any(iE) && max(max(abs(E(iE, :)))) <= 2*pi + 1e-3
    warning('Euler angles look like radians; converting to degrees.');
    E = E * 180 / pi;
end

%% Place every point on the grid from its coordinates
ix = round((X - min(X)) / XStep) + 1;
iy = round((Y - min(Y)) / YStep) + 1;
nx = max(ix); ny = max(iy);
lin = sub2ind([ny, nx], iy, ix);

phaseMap = zeros(ny, nx);
eulerMap = zeros(ny * nx, 3);
phaseMap(lin) = ph;
eulerMap(lin, :) = E;
missing = ny * nx - numel(unique(lin));
if missing > 0
    fprintf('  Note: %d grid positions have no data (treated as non-indexed).\n', missing);
end

indexed = ismember(phaseMap, hexPhase);

%% Orientation matrices and quaternions (vectorised, same formulas as
%  Gmatrix.m and Quaternion_Euler_Deg.m)
p1 = eulerMap(:, 1) * pi/180; P = eulerMap(:, 2) * pi/180; p2 = eulerMap(:, 3) * pi/180;
c1 = cos(p1); s1 = sin(p1); c2 = cos(p2); s2 = sin(p2); cP_ = cos(P); sP = sin(P);
g11 = c1.*c2 - s1.*s2.*cP_;  g12 = s1.*c2 + c1.*s2.*cP_;  g13 = s2.*sP;
g21 = -c1.*s2 - s1.*c2.*cP_; g22 = -s1.*s2 + c1.*c2.*cP_; g23 = c2.*sP;
g31 = s1.*sP;                g32 = -c1.*sP;               g33 = cP_;
G = [g11 g21 g31 g12 g22 g32 g13 g23 g33];   % column-major (= reshape(g,1,9))

c = cos(P/2); s = sin(P/2); sigma = 0.5*(p1 + p2); delta = 0.5*(p1 - p2);
Q = [c.*cos(sigma), s.*cos(delta), s.*sin(delta), c.*sin(sigma)];
Q(~indexed(:), :) = 0;

%% Output
ebsd.name     = name;
ebsd.file     = filename;
ebsd.ny = ny; ebsd.nx = nx;
ebsd.step = XStep; ebsd.stepY = YStep;
ebsd.phase = phaseMap;
ebsd.euler = reshape(eulerMap, ny, nx, 3);
ebsd.indexed = indexed;
ebsd.G = G;
ebsd.Q = Q;
ebsd.phases = phases;
ebsd.hexPhase = hexPhase;
ebsd.header = lines;

fprintf('  %s: %d x %d pixels, step %.4g um\n', name, nx, ny, XStep);
for k = 1:numel(phases)
    tag = '';
    if ismember(k, hexPhase), tag = '  <- used for macrozones'; end
    fprintf('    Phase %d: %-16s (%5.1f%% of map)%s\n', k, phases(k).name, ...
        100 * mean(phaseMap(:) == k), tag);
end
fprintf('    Non-indexed:  %5.1f%%\n', 100 * mean(phaseMap(:) == 0));
end

%% ---------------------------------------------------------------------
function v = mtr_readValue(lines, key)
v = NaN;
for k = 1:numel(lines)
    t = regexp(lines{k}, ['^' key '\s+([-+0-9.eE]+)'], 'tokens', 'once');
    if ~isempty(t), v = str2double(t{1}); return; end
end
end

function d = mtr_minStep(v)
u = unique(round(v * 1e4) / 1e4);
d = min(diff(u));
if isempty(d), d = 1; end
end

function R = mtr_DetectMacrozones(ebsd, P)
%MACROZONE_TOOL  Macrozone identification algorithm of
%   B. Fernandez Silva, M. Jackson, K. Fox, B.P. Wynne,
%   J. Appl. Cryst. 56 (2023) 737-749 (UoS_Macrozone_Tool.m), adapted to
%   work directly on the output of mtr_ReadCTF.
%
%   The logic is the same as in the original script (steps 1-7). Changes:
%   - The disorientation criterion (Trace_Disorientation) is evaluated for
%     the whole window at once with matrix products -> same result, much
%     faster on large maps.
%   - The quaternion sum of each macrozone is accumulated instead of being
%     recomputed from scratch at every step (same value).
%   - The 'macrocell' array of (max size x no. of macrozones), which runs out
%     of memory on large maps, is avoided.
%
%   P: struct with fields gridsize, critDisor (deg), critFraction (0-1),
%      areaThreshold (um^2)

if nargin < 2
    disp('mtr_DetectMacrozones needs data read with mtr_ReadCTF.');
    R = [];
    return
end

ny = ebsd.ny; nx = ebsd.nx;
Indexed = ebsd.indexed;
Stepsize = ebsd.step;
G = ebsd.G;  Qd = ebsd.Q;

gridsize = P.gridsize;
if mod(gridsize, 2) == 0
    error('gridsize must be ODD (the window is centred on the pixel).');
end
CDrad = P.critDisor * pi/180;
CF = P.critFraction;
thrT = 1 + 2*cos(CDrad);   % acos(0.5*(T-1)) <= CD  <=>  T >= 1 + 2cos(CD)
tolT = 1e-9;                % margin for re-checking borderline cases

% 12 hexagonal symmetry operators (same as SHex_Cell.mat)
r3 = sqrt(3)/2;
SH = {[1 0 0;0 1 0;0 0 1], [0.5 r3 0;-r3 0.5 0;0 0 1], [-0.5 r3 0;-r3 -0.5 0;0 0 1], ...
      [-1 0 0;0 -1 0;0 0 1], [-0.5 -r3 0;r3 -0.5 0;0 0 1], [0.5 -r3 0;r3 0.5 0;0 0 1], ...
      [1 0 0;0 -1 0;0 0 -1], [0.5 r3 0;r3 -0.5 0;0 0 -1], [-0.5 r3 0;r3 0.5 0;0 0 -1], ...
      [-1 0 0;0 1 0;0 0 -1], [-0.5 -r3 0;-r3 0.5 0;0 0 -1], [0.5 -r3 0;-r3 -0.5 0;0 0 -1]};
Sstack = vertcat(SH{:});                         % 36 x 3
symM = @(ga) reshape(permute(reshape(Sstack*ga, 3, 12, 3), [1 3 2]), 9, 12);
% T(:,s) = trace(gb * (S_s*ga)')  for every gb in the window

%% Size constants (identical to the original)
Size1 = round(gridsize/2);
Size2 = ny - round(gridsize/2);
Size3 = nx - round(gridsize/2);
Size4 = floor(gridsize/2);
Size5 = gridsize^2 - 1;
Size6 = ny - floor(gridsize/2);
Size7 = nx - floor(gridsize/2);

% linear offsets of the window
% Step 1 (order irrelevant) and Step 2 (a outer, b inner, as the original)
aa = kron(-Size4:Size4, ones(1, gridsize));   % row
bb = repmat(-Size4:Size4, 1, gridsize);       % column
offWin = aa(:) + bb(:)*ny;

tStart = tic;
%% STEP 1 - ismacro
ismacro = false(ny, nx);
for j = Size1:Size2
    for i = Size1:Size3
        c = j + (i-1)*ny;
        if Indexed(c) && ~ismacro(c)
            w = c + offWin;
            w = w(Indexed(w));
            ga = reshape(G(c, :), 3, 3);
            T = G(w, :) * symM(ga);
            ok = any(T >= thrT, 2);
            dud = any(abs(T - thrT) < tolT, 2);       % borderline cases: exact formula
            for q = find(dud)'
                ok(q) = mtr_exactTrace(reshape(G(w(q), :), 3, 3), ga, SH, CDrad);
            end
            Total = sum(ok);
            if (Total - 1)/Size5 >= CF
                ismacro(w(ok)) = true;
            end
        end
    end
end
fprintf('  Step 1 (MTR) completed: %.2f min\n', toc(tStart)/60);

%% STEP 2 - MacroNumber (region growing with mean orientation)
MacroNumber = zeros(ny, nx);
Mean_Orientation_Macrozones = zeros(0, 4);
Macronum = 0;
cap = 10000;
for j = Size1:Size2
    for i = Size1:Size3
        if ismacro(j, i) && MacroNumber(j, i) == 0
            Macronum = Macronum + 1;
            zone = zeros(cap, 2);  zone(1, :) = [j, i];
            dummy2 = 1; dummy3 = 1;
            AVG = zeros(3, 3); hasMean = false;
            qsum = zeros(1, 4); nq = 0;
            while true
                r0 = zone(dummy2, 1); c0 = zone(dummy2, 2);
                if r0 > Size4 && r0 < Size6 && c0 > Size4 && c0 < Size7
                    if hasMean
                        gc = AVG;
                    else
                        gc = reshape(G(r0 + (c0-1)*ny, :), 3, 3);
                    end
                    w = r0 + (c0-1)*ny + offWin;
                    cand = MacroNumber(w) == 0 & ismacro(w);
                    if any(cand)
                        wc = w(cand);
                        T = G(wc, :) * symM(gc);
                        okM = T >= thrT;
                        ok = any(okM, 2);
                        % trace of the FIRST symmetry that passes (as Trace_forQmean);
                        % the original leaves points with Theta = 0 out of the mean
                        [~, s1] = max(okM, [], 2);
                        Tfirst = T(sub2ind(size(T), (1:numel(wc))', s1));
                        valid = Tfirst < 3;
                        dud = any(abs(T - thrT) < tolT, 2) | (ok & abs(Tfirst - 3) < tolT);
                        for q = find(dud)'
                            [ok(q), valid(q)] = mtr_exactTrace(reshape(G(wc(q), :), 3, 3), gc, SH, CDrad);
                        end
                        if any(ok)
                            wc = wc(ok); valid = valid(ok);
                            MacroNumber(wc) = Macronum;
                            n = numel(wc);
                            if dummy3 + n > size(zone, 1)
                                zone = [zone; zeros(max(cap, n), 2)]; %#ok<AGROW>
                            end
                            [rr, cc] = ind2sub([ny nx], wc);
                            zone(dummy3+1:dummy3+n, :) = [rr, cc];
                            dummy3 = dummy3 + n;
                            if any(valid)
                                qsum = qsum + sum(Qd(wc(valid), :), 1);
                                nq = nq + sum(valid);
                            end
                        end
                    end
                    if nq > 0
                        Qm = qsum / norm(qsum);
                        AVG = mz_GMatrix_Quaternion(Qm(1), Qm(2), Qm(3), Qm(4));
                        hasMean = true;
                    end
                end
                if dummy2 >= dummy3
                    EULER = mz_EulerAngles_Gmatrix(AVG);
                    Mean_Orientation_Macrozones(Macronum, :) = [Macronum, EULER];
                    break
                end
                dummy2 = dummy2 + 1;
            end
        end
    end
end
fprintf('  Step 2 (macrozone numbering) completed: %.2f min - %d macrozones\n', ...
    toc(tStart)/60, Macronum);

%% STEP 3 - pixels of each macrozone (inside the evaluated region)
sub = MacroNumber(Size1:Size2, Size1:Size3);
[JJ, II] = find(sub > 0);
lab = sub(sub > 0);
JJ = JJ + Size1 - 1; II = II + Size1 - 1;
[~, ord] = sortrows([lab, JJ, II]);
JJ = JJ(ord); II = II(ord); lab = lab(ord);
nummac = accumarray(lab, 1, [Macronum 1]);
pix = mat2cell([JJ II], nummac, 2);

%% STEP 4 - colours
macrocell2 = cell(Macronum, 18);
colours = rand(Macronum, 3);
MacroNumberColor = zeros(ny, nx, 3);
regionMask = false(ny, nx); regionMask(Size1:Size2, Size1:Size3) = true;
m = MacroNumber > 0 & regionMask;
for k = 1:3
    tmp = zeros(ny, nx);
    tmp(m) = colours(MacroNumber(m), k);
    MacroNumberColor(:, :, k) = tmp;
end
for k = 1:Macronum
    macrocell2{k, 1} = k;
    macrocell2{k, 2} = pix{k};
    macrocell2{k, 12} = false;
end

%% STEP 5 - perimeters (convhull, boundary 0.5 and boundary 1)
for k = 1:Macronum
    XY = pix{k};
    if size(XY, 1) > 10
        X = XY(:, 1); Y = XY(:, 2);
        try
            [kk, vk] = convhull(X, Y);
            [j1, vj] = boundary(X, Y);
            [j2, vj2] = boundary(X, Y, 1);
            macrocell2{k, 3} = [X(kk), Y(kk)];  macrocell2{k, 4} = vk;
            macrocell2{k, 6} = [X(j1), Y(j1)];  macrocell2{k, 7} = vj;
            macrocell2{k, 9} = [X(j2), Y(j2)];  macrocell2{k, 10} = vj2;
            macrocell2{k, 12} = true;
        catch
        end
    end
end

%% STEP 6 - ellipses
for k = 1:Macronum
    try
        b1 = macrocell2{k, 6};
        if size(b1, 1) > 5
            a = mz_fitellipse(b1(:, 2), b1(:, 1)); a(5) = a(5) + pi/2; macrocell2{k, 8} = a;
        end
        ch = macrocell2{k, 3};
        if size(ch, 1) > 5
            a = mz_fitellipse(ch(:, 2), ch(:, 1)); a(5) = a(5) + pi/2; macrocell2{k, 5} = a;
        end
        b2 = macrocell2{k, 9};
        if size(b2, 1) > 5
            a = mz_fitellipse(b2(:, 2), b2(:, 1)); a(5) = a(5) + pi/2;
            macrocell2{k, 11} = a; macrocell2{k, 17} = a(2); macrocell2{k, 18} = a(1);
        end
    catch
    end
end
fprintf('  Steps 3-6 (perimeters and ellipses) completed: %.2f min\n', toc(tStart)/60);

%% STEP 7 - area threshold
Area_threshold = zeros(Macronum, 1);
for k = 1:Macronum
    if macrocell2{k, 12}
        Area_threshold(k) = macrocell2{k, 7} * Stepsize * Stepsize;
    end
end
removed = find(Area_threshold < P.areaThreshold);
MacroNumber_count = MacroNumber;
MacroNumber_count(ismember(MacroNumber_count, removed)) = 0;
ismacro_count = MacroNumber_count > 0;
MacroNumberColor_count = MacroNumberColor .* repmat(double(MacroNumber_count > 0), [1 1 3]);
kept = setdiff(1:Macronum, removed);
kept = kept(cellfun(@(v) isequal(v, true), macrocell2(kept, 12)));
kept = kept(:);

%% Statistics
nK = numel(kept);
ID = kept;
Area = zeros(nK, 1); Density = zeros(nK, 1); LongAx = zeros(nK, 1); ShortAx = zeros(nK, 1);
EllipseAngle = zeros(nK, 1); CentreX = zeros(nK, 1); CentreY = zeros(nK, 1);
okE = true(nK, 1);
for t = 1:nK
    k = kept(t);
    Area(t) = macrocell2{k, 7} * Stepsize^2;
    Density(t) = nummac(k) / macrocell2{k, 7} * 100;
    a = macrocell2{k, 11};
    if isempty(a) || numel(a) < 5 || any(~isfinite(a))
        okE(t) = false; continue
    end
    LongAx(t) = max(a(3), a(4)); ShortAx(t) = min(a(3), a(4));
    EllipseAngle(t) = a(5) * 180/pi; CentreX(t) = a(1); CentreY(t) = a(2);
end
Eu = Mean_Orientation_Macrozones(kept, 2:4);
if isempty(Eu), Eu = zeros(0, 3); end
DataTable = table(ID, nummac(kept), Area, sqrt(Area), sqrt(4*Area/pi), Density, ...
    ShortAx*Stepsize*2, LongAx*Stepsize*2, EllipseAngle, CentreX*Stepsize, CentreY*Stepsize, ...
    Eu(:, 1), Eu(:, 2), Eu(:, 3), ...
    'VariableNames', {'Macrozone', 'N_pixels', 'Area_um2', 'EqDiam_square_um', ...
    'EqDiam_circle_um', 'Density_pct', 'Short_axis_um', 'Long_axis_um', ...
    'Ellipse_angle_deg', 'Centre_X_um', 'Centre_Y_um', 'phi1_mean', 'PHI_mean', 'phi2_mean'});

MTR = 100 * sum(ismacro_count(:)) / (ny*nx);
MTR_alpha = 100 * sum(ismacro_count(:)) / max(1, sum(Indexed(:)));
De = sqrt(4*Area/pi); Dq = sqrt(Area); L = LongAx(okE)*Stepsize*2; Dn = Density;
est = @(v) [mean(v); max([v; NaN]); min([v; NaN]); std(v)];
values = [est(Area); est(Dq); est(De); est(Dn); est(L); MTR; MTR_alpha; nK; Macronum];
rowNames = {'Area_mean_um2', 'Area_max_um2', 'Area_min_um2', 'Area_std_um2', ...
    'EqDiam_square_mean_um', 'EqDiam_square_max_um', 'EqDiam_square_min_um', 'EqDiam_square_std_um', ...
    'EqDiam_circle_mean_um', 'EqDiam_circle_max_um', 'EqDiam_circle_min_um', 'EqDiam_circle_std_um', ...
    'Density_mean_pct', 'Density_max_pct', 'Density_min_pct', 'Density_std_pct', ...
    'LongAxis_mean_um', 'LongAxis_max_um', 'LongAxis_min_um', 'LongAxis_std_um', ...
    'MTR_pct_whole_map', 'MTR_pct_of_alpha_phase', 'N_macrozones_above_threshold', 'N_macrozones_total'};
Statistics = table(values, 'VariableNames', {'Value'}, 'RowNames', rowNames);

fprintf('  Completed: %.2f min - %d macrozones above threshold, MTR = %.1f%%\n', ...
    toc(tStart)/60, nK, MTR);

%% Result
R.P = P;
R.ismacro = ismacro;
R.MacroNumber = MacroNumber;
R.MacroNumberColor = MacroNumberColor;
R.MacroNumber_count = MacroNumber_count;
R.MacroNumberColor_count = MacroNumberColor_count;
R.ismacro_count = ismacro_count;
R.macrocell2 = macrocell2;
R.Macronum = Macronum;
R.nummac = nummac;
R.Area_threshold = Area_threshold;
R.Mean_Orientation_Macrozones = Mean_Orientation_Macrozones;
R.kept = kept;
R.DataTable = DataTable;
R.Statistics = Statistics;
R.MTR = MTR;
R.Fraction_Macro = sum(ismacro(:)) / ((Size2-Size1)*(Size3-Size1));
end

%% -------------------------------------------------------------------------
function [got, valid] = mtr_exactTrace(gb, ga, SH, disvalue)
% Exact replica of Trace_Disorientation / Trace_forQmean (for numerical
% ties): got = criterion met; valid = real Theta different from 0.
got = false; valid = false;
for s = 1:12
    T = trace(gb*(SH{s}*(ga)).');
    if (acos(0.5*(T-1))) <= disvalue
        got = true;
        valid = real(acosd(0.5*(T-1))) ~= 0;
        return
    end
end
end

%% ===== Local copies of the 'Stored Functions' helpers (renamed mz_* so they
%  never shadow functions of MTEX or of other scripts on the MATLAB path) =====

function [G] = mz_GMatrix_Quaternion(q0,q1,q2,q3)
qq = q0^2-(q1^2+q2^2+q3^2);
G = zeros (3,3);
G (1,1) = qq+2*(q1^2);
G (1,2) = 2*(q1*q2-q0*q3);
G (1,3) = 2*(q1*q3+q0*q2);
G (2,1) = 2*(q1*q2+q0*q3);
G (2,2) = qq+2*(q2^2);
G (2,3) = 2*(q2*q3-q0*q1);
G (3,1) = 2*(q1*q3-q0*q2);
G (3,2) = 2*(q2*q3+q0*q1);
G (3,3) = qq+2*(q3^2);
 G = G.';  
end

function [ EA ] = mz_EulerAngles_Gmatrix( G )
 C= 1/(1-(G(3,3))^2)^1/2;
    EA = zeros(1,3);
    if  (abs(G(3,3))==1)
     EA (1,1)= atan2(G(1,2),G(1,1));
     EA (1,2)= acos (G(3,3));
     EA (1,3)=0;
    else 
     EA (1,1)= atan2((G(3,1)*C),-(G(3,2)*C));
     EA (1,2)= acos (G(3,3));
     EA (1,3)= atan2((G(1,3)*C),(G(2,3)*C));
    end
    if EA(1,1)<0 % Corrections
        EA(1,1)=EA(1,1)+2*pi;
    end
    if EA(1,3)<0
        EA(1,3)=EA(1,3)+2*pi;
    end
    EA=EA*180/pi;
end

function a = mz_fitellipse(X,Y)
mx = mean(X);
my = mean(Y);
sx = (max(X)-min(X))/2;
sy = (max(Y)-min(Y))/2; 
x = (X-mx)/sx;
y = (Y-my)/sy;
 x = x(:);
 y = y(:);
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];
S = D'*D;
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
if 0
 [gevec, geval] = eig(S,C);
 I = find(real(diag(geval)) < 1e-8 & ~isinf(diag(geval)));
 A = real(gevec(:,I));
else
  tmpA = S(1:3,1:3); 
  tmpB = S(1:3,4:6); 
  tmpC = S(4:6,4:6); 
  tmpD = C(1:3,1:3);
  tmpE = inv(tmpC)*tmpB';
  [evec_x, eval_x] = eig(inv(tmpD) * (tmpA - tmpB*tmpE));
  I = find(real(diag(eval_x)) < 1e-8 & ~isinf(diag(eval_x)));
  A = real(evec_x(:,I));
  evec_y = -tmpE * A;
  A = [A; evec_y];
par = [
  A(1)*sy*sy,   ...
      A(2)*sx*sy,   ...
      A(3)*sx*sx,   ...
      -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
      -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
      A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
      - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
      + A(6)*sx*sx*sy*sy   ...
      ]';
thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;
Ao = par(6);
Au =   par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;
tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;
uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;
Ru = -wCentre./Auu;
Rv = -wCentre./Avv;
Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);
a = [uCentre, vCentre, Ru, Rv, thetarad,];
end
end
