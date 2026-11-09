%% Ti64_Schmid_LusterMorris_MTEX61.m
% SCHMID FACTOR + LUSTER-MORRIS m' FOR ALPHA Ti-6Al-4V (2D EBSD / MTEX 6.x)
%
% This script:
%  (1) imports an Oxford .ctf EBSD map and selects the HCP alpha phase;
%  (2) reconstructs alpha grains, retaining grain-boundary misorientations;
%  (3) maps maximum |Schmid factor| for basal, prismatic <a>, pyramidal <c+a>;
%  (4) maps Luster-Morris m' on alpha/alpha grain boundaries:
%        a. load-selected basal -> load-selected basal
%        b. load-selected prismatic -> load-selected prismatic
%        c. geometric BEST basal -> basal (all candidate pairs)
%        d. geometric BEST prismatic -> prismatic (all candidate pairs)
%  (5) plots m' vs. misorientation, and exports figures and CSV files.
%
% IMPORTANT SCIENTIFIC LIMITATIONS
% - The load direction is expressed in the EBSD SAMPLE reference frame.
% - Schmid factors are geometric, not actual local stresses or slip activity.
% - m' indicates geometric slip-transfer compatibility, NOT observed transmission.
% - This does NOT simulate crystal plasticity, load shedding, dwell fatigue,
%   grain-level stresses, or phase alpha/beta slip transfer.
% - 2D EBSD produces 2D grain-boundary traces, NOT 3D boundary plane normals.
% - Parameters for grain segmentation MUST be checked against map resolution.
%
% MTEX background:
% https://mtex-toolbox.github.io/SlipTransmission.html
% https://mtex-toolbox.github.io/SchmidFactor.html
% https://mtex-toolbox.github.io/SlipSystems.html
%
% Authoring convention: m' = |(n1.n2)*(b1.b2)| in a COMMON sample frame.
% Requirements: MATLAB R2019b+ and MTEX 6.x/compatible releases on path.

%% ======================== EDIT OPTIONS HERE =============================
opt.filename      = '';       % '' -> ask which .ctf file to analyze
opt.phaseName     = '';       % '' -> detect alpha automatically or ask
opt.loadXYZ       = [-1 0 0]; % assumed UNAXIAL loading axis in EBSD frame
opt.segAngleDeg   = 5;        % grain boundary minimum misorientation (degrees)
opt.minPixels     = 5;        % remove smaller indexed clusters
opt.smoothSteps   = 3;        % visual smoothing, 0 means none
opt.lineWidthMP   = 2.0;      % grain boundary line width in m' maps
opt.showFigures   = true;     % make visible figures (not merely save them)
opt.outputDir     = '';       % '' -> subfolder alongside input CTF
%% =======================================================================

if exist('EBSD','class') ~= 8 || exist('slipSystem','class') ~= 8
    error('MTEX classes not found. Run your MTEX startup script first.');
end

if isempty(opt.filename)
    [shortName, sourceDir] = uigetfile({'*.ctf','Oxford EBSD (*.ctf)'}, ...
        'Select Ti-6Al-4V EBSD map');
    if isequal(shortName,0)
        error('No CTF file selected. Analysis cancelled.');
    end
    fileName = fullfile(sourceDir,shortName);
else
    fileName = opt.filename;
    [sourceDir,~,~] = fileparts(fileName);
    if isempty(sourceDir), sourceDir = pwd; end
end
if exist(fileName,'file') ~= 2
    error('CTF file not found: %s',fileName);
end
if isempty(opt.outputDir)
    outDir = fullfile(sourceDir,'Ti64_Schmid_LusterMorris_Output');
else
    outDir = opt.outputDir;
end
if ~exist(outDir,'dir'), mkdir(outDir); end

L = vector3d(opt.loadXYZ(1),opt.loadXYZ(2),opt.loadXYZ(3));
if norm(opt.loadXYZ) < eps, error('opt.loadXYZ cannot be zero.'); end
L = L ./ norm(L);

fprintf('\n===== Ti64 SCHMID + LUSTER-MORRIS =====\n');
fprintf('Input: %s\n',fileName);
fprintf('Assumed loading axis (sample coordinates): [%g %g %g]\n',opt.loadXYZ);
fprintf('Note: reversing the loading AXIS does not change |Schmid| or m''.\n');

%% 1) Import CTF using a consistent spatial/Euler reference frame
% IMPORTANT: if your existing research code uses a different import alignment,
% adopt THAT SAME convention here before comparing the maps with earlier work.
ebAll = EBSD.load(fileName,'interface','ctf', ...
    'convertEuler2SpatialReferenceFrame');
disp(ebAll);

% Determine the name of alpha phase; do not accidentally combine alpha/beta.
allPhaseNames = ebAll.mineralList;
if ischar(allPhaseNames), allPhaseNames = cellstr(allPhaseNames); end
if isstring(allPhaseNames), allPhaseNames = cellstr(allPhaseNames); end
valid = ~strcmpi(allPhaseNames,'notIndexed') & ~cellfun(@isempty,allPhaseNames);
phaseNames = allPhaseNames(valid);

if isempty(opt.phaseName)
    namesLower = lower(string(phaseNames));
    isAlphaName = contains(namesLower,'alpha') | contains(namesLower,'hcp') | ...
        contains(namesLower,char(945));
    if nnz(isAlphaName)==1
        phaseName = phaseNames{find(isAlphaName,1)};
    elseif numel(phaseNames)==1
        phaseName = phaseNames{1};
        warning('Only one indexed phase present. Verify this is HCP alpha!');
    else
        fprintf('\nAvailable indexed phase names:\n');
        for i = 1:numel(phaseNames)
            fprintf(' %d) %s\n',i,phaseNames{i});
        end
        k = input('Enter the number for HCP alpha phase: ');
        if isempty(k) || ~isscalar(k) || k<1 || k>numel(phaseNames) || k~=round(k)
            error('Invalid alpha-phase selection. Set opt.phaseName manually.');
        end
        phaseName = phaseNames{k};
    end
else
    phaseName = opt.phaseName;
end

ebAlpha = ebAll(phaseName);
if isempty(ebAlpha)
    error('No EBSD measurements found for phase: %s',phaseName);
end
cs = ebAlpha.CS;
fprintf('Selected phase: %s\n',phaseName);
disp(cs);

%% 2) Grain reconstruction from indexed alpha pixels only
[grains,ebAlpha] = calcGrains(ebAlpha,'angle', ...
    opt.segAngleDeg*degree,'minPixel',opt.minPixels);
if isempty(grains), error('No alpha grains after segmentation.'); end
if opt.smoothSteps > 0
    % For orientation-based boundary comparisons, labels remain unchanged.
    grains = smoothBoundary(grains,opt.smoothSteps);
end

gB = grains.boundary('indexed');
if isempty(gB)
    error('No alpha/alpha boundaries. Check region, resolution, segmentation.');
end

% gB.grainId consists of GRAIN IDs, not necessarily indices into grains.
segmentPairIDs = sort(gB.grainId,2);
[pairIDs,~,pairRowForSegment] = unique(segmentPairIDs,'rows');
ind1 = grains.id2ind(pairIDs(:,1));
ind2 = grains.id2ind(pairIDs(:,2));
if any(~isfinite(ind1)) || any(~isfinite(ind2)) || ...
        any(ind1<1) || any(ind2<1)
    error('Some grain-boundary IDs are not present in the grains object.');
end
fprintf('Alpha grains: %d; boundary segments: %d; unique grain pairs: %d\n', ...
    length(grains),length(gB),size(pairIDs,1));

%% 3) Slip-system families, with both signs treated as geometrically equivalent
% NOTE: do not rank different families by Schmid factor alone; their CRSS
% differ. Results for each family are kept SEPARATE here.
sBas = slipSystem.basal(cs).symmetrise('antipodal');
sPri = slipSystem.prismaticA(cs).symmetrise('antipodal');
sPyr = slipSystem.pyramidalCA(cs).symmetrise('antipodal');

oriG = grains.meanOrientation;
Lcrystal = inv(oriG) * L;
[SFbas, iBas] = max(abs(sBas.SchmidFactor(Lcrystal)),[],2);
[SFpri, iPri] = max(abs(sPri.SchmidFactor(Lcrystal)),[],2);
[SFpyr, iPyr] = max(abs(sPyr.SchmidFactor(Lcrystal)),[],2);

% The selected systems refer to max absolute Schmid for EACH family.
sBasSelected = oriG .* sBas(iBas);
sPriSelected = oriG .* sPri(iPri);

%% 4) Luster-Morris m' at unique alpha/alpha boundaries
% (i) load-selected systems: Schmid-selected on each side, within a family.
mPselectedBas = abs(mPrime(sBasSelected(ind1),sBasSelected(ind2)));
mPselectedPri = abs(mPrime(sPriSelected(ind1),sPriSelected(ind2)));

% (ii) BEST geometrically compatible pairing (no load-selection condition).
% Each i,j pair of basal/prismatic systems is checked independently.
mPbestBas = bestPairwiseMPrime(oriG(ind1),oriG(ind2),sBas);
mPbestPri = bestPairwiseMPrime(oriG(ind1),oriG(ind2),sPri);

mPselectedBas = min(1,max(0,mPselectedBas(:)));
mPselectedPri = min(1,max(0,mPselectedPri(:)));
mPbestBas     = min(1,max(0,mPbestBas(:)));
mPbestPri     = min(1,max(0,mPbestPri(:)));

% Symmetry-reduced grain misorientation; not equal to slip-system mismatch.
misDeg = angle(oriG(ind1),oriG(ind2)) ./ degree;
misDeg = misDeg(:);

% Inclination between alpha c-axis / basal normal and loading AXIS.
cAxis = Miller(0,0,0,1,cs,'hkl');
thetaC = angle(oriG * cAxis,L) ./ degree;
thetaC = min(thetaC(:),180-thetaC(:));

%% 5) Grain-coloured Schmid maps (one independent colorbar per figure)
vis = 'on';
if ~opt.showFigures, vis = 'off'; end
h = plotGrainQuantity(grains,SFbas,'Schmid |m| - basal <a>', ...
    '|Schmid|, basal (0 to 0.5)',[0 0.5],vis);
writeFig(h,outDir,'01_Schmid_Basal');

h = plotGrainQuantity(grains,SFpri,'Schmid |m| - prismatic <a>', ...
    '|Schmid|, prismatic (0 to 0.5)',[0 0.5],vis);
writeFig(h,outDir,'02_Schmid_Prismatic');

h = plotGrainQuantity(grains,SFpyr,'Schmid |m| - pyramidal <c+a>', ...
    '|Schmid|, pyramidal c+a (0 to 0.5)',[0 0.5],vis);
writeFig(h,outDir,'03_Schmid_PyramidalCA');

%% 6) Boundary-coloured m' maps
h = plotBoundaryQuantity(grains,gB,mPselectedBas(pairRowForSegment), ...
    "m' - Schmid-selected basal / basal",'m'' basal selected',[0 1], ...
    opt.lineWidthMP,vis);
writeFig(h,outDir,'04_mPrime_Selected_Basal');

h = plotBoundaryQuantity(grains,gB,mPselectedPri(pairRowForSegment), ...
    "m' - Schmid-selected prismatic / prismatic",'m'' prismatic selected', ...
    [0 1],opt.lineWidthMP,vis);
writeFig(h,outDir,'05_mPrime_Selected_Prismatic');

h = plotBoundaryQuantity(grains,gB,mPbestBas(pairRowForSegment), ...
    "m' - maximum geometrical basal / basal",'max m'' basal / basal', ...
    [0 1],opt.lineWidthMP,vis);
writeFig(h,outDir,'06_mPrime_Best_Basal');

h = plotBoundaryQuantity(grains,gB,mPbestPri(pairRowForSegment), ...
    "m' - maximum geometrical prismatic / prismatic",'max m'' prism / prism', ...
    [0 1],opt.lineWidthMP,vis);
writeFig(h,outDir,'07_mPrime_Best_Prismatic');

%% 7) Histograms and correlation against misorientation
h = figure('Name','mPrime distributions','Color','w','Visible',vis);
histogram(mPselectedBas,0:0.05:1,'Normalization','probability', ...
    'DisplayStyle','stairs','LineWidth',1.8);
hold on
histogram(mPbestBas,0:0.05:1,'Normalization','probability', ...
    'DisplayStyle','stairs','LineWidth',1.8);
hold off
xlabel('Luster-Morris m'''); ylabel('Fraction of unique grain pairs');
title('Basal compatibility: load-selected vs maximum geometric');
legend({'Selected by Schmid','Maximum over basal pairs'},'Location','best');
grid on; xlim([0 1]);
writeFig(h,outDir,'08_mPrime_Histogram');

h = figure('Name','mPrime versus misorientation','Color','w','Visible',vis);
scatter(misDeg,mPbestBas,8,mPselectedBas,'filled');
xlabel('alpha/alpha misorientation (degrees)');
ylabel('Maximum basal m'' (geometric)');
title('Boundary pairs: colour = load-selected basal m''');
colormap(parula(256)); climCompat([0 1]);
cb = colorbar; ylabel(cb,'Selected basal m''');
ylim([0 1]); grid on;
writeFig(h,outDir,'09_mPrime_vs_Misorientation');

%% 8) Export numerical data; one boundary row = one UNIQUE grain pair
cg = grains.centroid;
grainTable = table(grains.id(:),cg.x(:),cg.y(:),grains.numPixel(:), ...
    thetaC,SFbas(:),SFpri(:),SFpyr(:), ...
    'VariableNames',{'GrainID','CentroidX','CentroidY','NumPixels', ...
    'Caxis_Load_Angle_deg','Schmid_Basal','Schmid_Prismatic', ...
    'Schmid_PyramidalCA'});
writetable(grainTable,fullfile(outDir,'Grain_Schmid.csv'));

boundaryTable = table(pairIDs(:,1),pairIDs(:,2),misDeg, ...
    mPselectedBas,mPselectedPri,mPbestBas,mPbestPri, ...
    'VariableNames',{'GrainID_1','GrainID_2','Misorientation_deg', ...
    'mPrime_Selected_Basal','mPrime_Selected_Prismatic', ...
    'mPrime_Best_Basal','mPrime_Best_Prismatic'});
writetable(boundaryTable,fullfile(outDir,'Boundary_mPrime.csv'));

fprintf('\n--- NUMERICAL RESULTS (unique grain pairs) ---\n');
fprintf('Schmid basal:   mean = %.3f, max = %.3f\n',mean(SFbas),max(SFbas));
fprintf('Schmid prism:   mean = %.3f, max = %.3f\n',mean(SFpri),max(SFpri));
fprintf('Schmid pyr CA:  mean = %.3f, max = %.3f\n',mean(SFpyr),max(SFpyr));
fprintf('m'' selected basal: median = %.3f\n',median(mPselectedBas));
fprintf('m'' best basal:     median = %.3f\n',median(mPbestBas));
fprintf('Results written to: %s\n',outDir);
fprintf('REMINDER: NO local stress redistribution or load shedding is calculated.\n');

%% ============================ LOCAL FUNCTIONS ==========================
function result = bestPairwiseMPrime(ori1,ori2,sS)
% Max m' across every slip-system pair of this family for each grain pair.
    nPairs = length(ori1);
    result = zeros(nPairs,1);
    for ia = 1:length(sS)
        ssA = ori1 .* sS(ia);
        for ib = 1:length(sS)
            ssB = ori2 .* sS(ib);
            candidate = abs(mPrime(ssA,ssB));
            result = max(result,candidate(:));
        end
    end
end

function h = plotGrainQuantity(grains,values,figTitle,barTitle,rng,visibility)
    h = figure('Color','w','Name',figTitle,'Visible',visibility);
    plot(grains,values(:),'figSize','large');
    hold on;
    plot(grains.boundary,'lineWidth',0.35,'lineColor',[0.18 0.18 0.18]);
    hold off;
    colormap(parula(256));
    climCompat(rng);
    mtexColorbar('title',barTitle);
    title(figTitle,'Interpreter','none');
    axis equal tight;
    drawnow;
end

function h = plotBoundaryQuantity(grains,gB,segmentValues,figTitle,barTitle,rng,lw,visibility)
    h = figure('Color','w','Name',char(figTitle),'Visible',visibility);
    plot(grains,'FaceColor',[0.88 0.88 0.88],'figSize','large');
    hold on;
    plot(gB,segmentValues(:),'lineWidth',lw);
    hold off;
    colormap(parula(256));
    climCompat(rng);
    mtexColorbar('title',barTitle);
    title(char(figTitle),'Interpreter','none');
    axis equal tight;
    drawnow;
end

function climCompat(limits)
% caxis in MTEX 6.x; clim in newer MATLAB. caxis is most portable here.
    caxis(limits);
end

function writeFig(h,d,name)
    drawnow;
    savefig(h,fullfile(d,[name '.fig']));
    try
        exportgraphics(h,fullfile(d,[name '.png']),'Resolution',250);
    catch ME
        warning('exportgraphics failed (%s); using print instead.',ME.message);
        print(h,fullfile(d,[name '.png']),'-dpng','-r250');
    end
end
