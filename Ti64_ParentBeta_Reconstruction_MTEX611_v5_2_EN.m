%% Ti-6Al-4V: prior-beta reconstruction, MTEX 6.1.1 (v5.2)
% v5.1 = v5.0 with the load arrow of the Schmid / stress maps (F01-F13) drawn in a
%        thin band between the title and the map instead of on top of the map.
%        + Schmid factors per MTR / macrozone, linked to Results.mat of
%          RUN_MACROZONES for the SAME .ctf (F00 extra rows, F00b table, F14-F15).
% v5.2 = MTR detection (Fernandez Silva et al., J. Appl. Cryst. 56 (2023) 737-749)
%        built in (SFopt.macrozoneResults = 'run'), MTR vs non-MTR Schmid
%        distributions (F16) and in-plane load-direction sweep per macrozone (F17).
% v5.0 = v4.9 with the figures de-duplicated and the PF kernel tied to Channel 5.
% NO change to the reconstruction, BOR, coverage, fragmentation or Schmid logic.
% NOT TESTED in MATLAB (edited without MATLAB access).
%
% FIGURES (file names and window titles use the category letter)
%  A  Measured map:        A1 alpha IPF-Z, A2 PHASE MAP (alpha / beta / zero solution)
%  B  Reconstruction:      B1/B2 BOR fit (initial/used), B3 beta before merge, B4 final
%                          beta IPF-Z, B5 beta boundaries on alpha, B6/B7 per-pixel beta + fit
%  C  Coverage/fragments:  C1 coverage map, C2 fragmentation, C3 beta ECD histograms
%  D  Textures (PF/IPF):   D1 WHOLE-MAP alpha PF (compare with Channel 5), D2 measured beta PF
%                          (only with showRedundantPF), D3 measured beta IPF X/Y/Z, D4 six PFs
%                          alpha vs reconstructed beta, D5 nine PFs, D6/D7 reconstructed beta PF
%                          (pixel / grain-area weighted), D8 scatter IPF-Z, D9/D10 repeats of D4
%                          (only with showRedundantPF), D0 CSV halfwidth check
%  E  BOR validation:      E1-E4 (+ E0 CSVs)
%  F  Schmid factor:       F01-F13 maps/histograms (+ F00 CSVs)
%  G  Convergence:         G1 (optional)
%
% PF KERNEL AND CHANNEL 5 (measured on SAMPLE015, 134025 alpha px)
%  MTEX halfwidth h reproduces Channel 5 "half width" 2h  (MTEX 5 deg = Channel 5
%  10 deg; MTEX 10 deg = Channel 5 20 deg). Both use ~Gaussian bell kernels, but
%  Channel 5 quotes ~2x the MTEX value. Set only pfHalfwidth (section 1).
%  Descriptive equivalence, not Oxford documentation. Channel 5 shows the RAW CTF
%  frame; MTEX here uses convertEuler2SpatialReferenceFrame, so compare MAXIMA,
%  not peak azimuths, unless you rotate 180 deg about X.
%
% A reconstructed beta orientation is an INFERENCE from the BOR; it does not
% prove that each alpha grain came from a single recoverable beta parent.

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
alphaPhase = 'Ti-Hex';
betaPhase  = 'Titanium cubic';

% --- Main-run parameters (a starting point, NOT "the truth") ---------------
P = struct();
P.childThreshold     = 1.5*degree; % alpha child-grain segmentation
P.minPixels          = 3;          % remove alpha grains < minPixels px, re-segment
P.graphThreshold     = 2.5*degree; % calcVariantGraph 'threshold'
P.graphTolerance     = 2.5*degree; % calcVariantGraph 'tolerance' ([] = not passed)
P.graphIterations    = 3;
P.minProb            = 0.50;
P.mergeThreshold     = 5*degree;
P.maxInclusionPixels = 10;
P.refineBOR          = false;      % true: fit the experimental BOR (calcParent2Child)
P.refQuantile        = 0.60;

% --- Your ORIGINAL parameters, for direct comparison -----------------------
Pbase = P;
Pbase.childThreshold  = 1.5*degree;
Pbase.minPixels       = 1;
Pbase.graphThreshold  = 1.5*degree;
Pbase.graphTolerance  = [];        % your code did not pass 'tolerance'
Pbase.minProb         = 0.50;
Pbase.mergeThreshold  = 5*degree;

% --- Region and resolution -------------------------------------------------
roiRect          = [];  % [x0 y0 width height] in CTF units (um); [] = whole map
downsampleFactor = 1;   % >1: subsample the grid (every k-th pixel)

% --- Optional modules (expensive) -----------------------------------------
runBaselineComparison = true;
runSensitivity        = false;
sensMode              = 'oneAtATime';   % 'oneAtATime' or 'fullGrid' (108 runs)
runConvergence        = false;
convFactors           = [1 2 3 4 6 8];  % subsampling factors of the same map

% --- Heuristic reliability criteria (ADJUSTABLE, not a standard) -----------
Q = struct();
Q.minPxPerAlphaDiameter  = 5;    % pixels across a typical alpha diameter
Q.maxMedianBOR_deg       = 2.5;
Q.maxP90BOR_deg          = 5;
Q.minMedianFragmentPx    = 10;
Q.maxPctFragLE3          = 20;   % % of beta grains with <= 3 px
Q.minInteriorGrains      = 30;
Q.minFovOverBetaECD      = 5;    % shorter map side / P90 of ALL beta ECDs
Q.maxPctEdge             = 30;   % % of beta grains touching edge or holes
Q.minPlausibleBetaECD_um = 3;    % beta grains smaller than this are implausible
Q.maxPctImplausible      = 20;   % % of grains allowed below that size

%% PF kernel (MTEX) <-> Channel 5 "half width"
% Channel 5 half width = 2 x MTEX halfwidth (see header).  pfHalfwidth is the ONLY
% kernel parameter used by figures 10, 11, 12, 25, 26, 27, 28 and the CSVs.
ch5HalfwidthTarget_deg = 10;                      % value you set in Channel 5
pfHalfwidth = (ch5HalfwidthTarget_deg/2)*degree;  % -> 5 deg in MTEX
runHalfwidthCheck  = true;                  % PF maxima vs halfwidth (whole-map alpha)
halfwidthList_deg  = [2.5 5 7.5 10 15];     % MTEX halfwidths to scan
ch5Ref0002         = [];                    % e.g. 18.86 = Channel 5 {0002} max; [] = skip
showRedundantPF    = false;                 % true: also figs D2, D9, D10 (repeat D5 / D4)
pfAlphaPopulation  = 'whole';   % 'whole' = ALL indexed alpha (like fig D1, same as Channel 5); 'matched' = only px with reconstructed beta
runScatterPF      = true;                  % D8: scatter IPF-Z (alpha / beta / reconstructed beta)
scatterMaxPoints   = 15000;                 % max dots per panel (display only)

% --- Optional Schmid-factor analysis of the MEASURED alpha phase -----------
% This is a geometric uniaxial-loading scenario, NOT a stress/strain solver.
% The arrow points left on the EBSD map: -X in the specimen/map frame.
% Verify that the EBSD X axis actually aligns with the physical radial axis.
runSchmidAnalysis = true;
% New optional quantitative BOR check of matched alpha / reconstructed beta.
runBORPairwiseCheck = true;
BORopt = struct();
BORopt.chunkSize = 8000;            % limit MTEX memory consumption
BORopt.maxOverlayPairs = 2000;     % representative pairs displayed in Fig. 18
BORopt.maxPairLines = 60;          % short lines connecting corresponding poles
BORopt.exportPixelCSV = false;     % true: large CSV; MAT always saved
BORopt.randomSeed = 11;           % reproducible PF plot subsampling
BORopt.toleranceDeg = [2 5];      % BOR angular thresholds [deg]
SFopt = struct(); %% DWELL FATIGUE LOADING CONDITIONS

% Loading direction in EBSD reference frame
SFopt.loadAxisXYZ = [-1 0 0];

% Yield strength (illustrative value)
yieldStrengthMPa = 900;

% Maximum applied stress during dwell fatigue
dwellStressRatio = 0.95;

% Nominal stress during dwell
SFopt.nominalStressMPa = ...
    dwellStressRatio * yieldStrengthMPa;
SFopt.crssMPa = [NaN NaN NaN];   % unknown: do not predict slip activation
SFopt.exportPixelCSV = false; % true: detailed potentially very large CSV
SFopt.candidateMTRroiRect = []; % optional [x_min y_min width height], just a REGION, not automatic MTR detection
% Link to the macrozone (MTR) analysis of the SAME .ctf (RUN_MACROZONES):
%   'run'  = detect the MTRs HERE, on the same .ctf (built-in macrozone tool,
%            Fernandez Silva et al. 2023); saved as MTR_Results.mat
%   'auto' = use <sample>_Macrozones_*\Results.mat next to the .ctf (newest)
%   'ask'  = choose the Results.mat in a window
%   ''     = do not link (no MTR rows / tables)
%   or a full path to a Results.mat
SFopt.macrozoneResults = 'run';
% MTR detection parameters (used with 'run'; paper defaults)
SFopt.MTR.gridsize      = 'auto';  % odd number of pixels, or 'auto' = ~windowUm window
SFopt.MTR.windowUm      = 50;
SFopt.MTR.critDisor     = 20;      % deg
SFopt.MTR.critFraction  = 0.3;
SFopt.MTR.areaThreshold = 10000;   % um^2
% Sensitivity of each macrozone to the in-plane load direction (F17)
SFopt.loadSweep        = true;
SFopt.loadSweepStepDeg = 5;        % 0:step:180 deg from +X, in the map plane
SFopt.cAxisClassDeg = [30 60];  % c-axis-to-load angle: <30 'hard', 30-60 'intermediate', >60 'soft'
% Note: +/-X give identical |Schmid factor| under uniaxial loading.

%% 2. Input file
[fileName,filePath] = uigetfile({'*.ctf','Oxford / AZtec CTF (*.ctf)'}, ...
    'Select the Ti-6Al-4V CTF');
if isequal(fileName,0), error('No file selected.'); end
fname = fullfile(filePath,fileName);
[~,sampleName] = fileparts(fname);
outputFolder = fullfile(filePath,[sampleName '_PriorBeta_MTEX611_v4_BOR']);
if ~isfolder(outputFolder), mkdir(outputFolder); end
fprintf('\nINPUT:  %s\nOUTPUT: %s\n',fname,outputFolder);

%% 3. Load CTF, ROI and subsampling
ebsdRaw = EBSD.load(fname,'interface','ctf','convertEuler2SpatialReferenceFrame');
disp(ebsdRaw);

if ~isempty(roiRect)
    ebsdRaw = ebsdRaw(inpolygon(ebsdRaw,roiRect));
    fprintf('ROI applied: [%s] -> %d points\n',num2str(roiRect),length(ebsdRaw));
end
stepNative = estimateStep(ebsdRaw);
if downsampleFactor > 1
    ebsdRaw = downsampleEBSD(ebsdRaw,downsampleFactor);
end
stepUm = estimateStep(ebsdRaw);
try, scanUnit = ebsdRaw.scanUnit; catch, scanUnit = 'unknown'; end
fprintf('Spatial unit of the map: %s | native step = %.4f | step used = %.4f\n', ...
    scanUnit,stepNative,stepUm);
if downsampleFactor > 1 && abs(stepUm/stepNative - downsampleFactor) > 0.05*downsampleFactor
    warning('Step after subsampling (%.4f) is not %dx the native step (%.4f): check.', ...
        stepUm,downsampleFactor,stepNative);
end

%% 4. Phases and physical indexed area
try
    alphaData = ebsdRaw(alphaPhase);
    betaData  = ebsdRaw(betaPhase);
    CSalpha = alphaData.CS;  CSbeta = betaData.CS;
catch
    disp(ebsdRaw);
    error('Phase names not found: edit alphaPhase/betaPhase.');
end
if isempty(alphaData), error('No indexed alpha pixels.'); end
mapW = max(ebsdRaw.x) - min(ebsdRaw.x) + stepUm;
mapH = max(ebsdRaw.y) - min(ebsdRaw.y) + stepUm;
mapMin = min(mapW,mapH);
fprintf('Indexed alpha pixels: %d | measured beta: %d\n',length(alphaData),length(betaData));
fprintf('Indexed alpha area = %.5g mm^2 (assuming um units) | map %.1f x %.1f\n', ...
    length(alphaData)*stepUm^2*1e-6,mapW,mapH);
if ~isempty(betaData)
    warning(['%d MEASURED beta pixels exist. job.parentGrains includes them; the ', ...
        'fragment statistics of this script exclude them.'],length(betaData));
end

%% 5. Main reconstruction
R = runReconstruction(ebsdRaw,alphaPhase,betaPhase,P,stepUm);
if ~R.ok
    error('Reconstruction produced no result: %s',R.reason);
end
fprintf('\nAlpha grains assigned after voting: %.1f%%\n',R.alphaGrainsAssigned_pct);

%% 6. Reconstruction figures (kept from the original script)
alphaKey = ipfHSVKey(CSalpha);  alphaKey.inversePoleFigureDirection = zvector;
betaKey  = ipfHSVKey(CSbeta);   betaKey.inversePoleFigureDirection  = zvector;

%% 6a. Phase map (alpha / beta / zero solution)
try
    plotPhaseMap(ebsdRaw,alphaData,betaData,sampleName,stepUm,outputFolder);
catch MEpm
    warning('Phase map failed: %s',MEpm.message);
end

try
    figure('Name','A1 Alpha IPF-Z','Color','white');
    plot(alphaData,alphaKey.orientation2color(alphaData.orientations));
    set(gca,'Color',[0.12 0.12 0.12]);
    title([sampleName ' | measured alpha IPF-Z'],'Interpreter','none');
    saveFigure(outputFolder,'A1_Alpha_IPFz.png');

    if ~isempty(R.gbFit0)
        figure('Name','B1 Initial BOR fit','Color','white'); histogram(R.gbFit0,50);
        xlabel('BOR misfit per boundary [deg]'); ylabel('Neighbour pairs');
        title('Ideal BOR: fit on alpha/alpha boundaries');
        saveFigure(outputFolder,'B1_BORfit_Initial.png');
    end
    if P.refineBOR && ~isempty(R.gbFit)
        figure('Name','B2 BOR fit used','Color','white'); histogram(R.gbFit,50);
        xlabel('BOR misfit per boundary [deg]'); ylabel('Neighbour pairs');
        title('BOR used in the reconstruction');
        saveFigure(outputFolder,'B2_BORfit_Used.png');
    end

    figure('Name','B3 Beta before merging','Color','white');
    plot(R.parentBefore,betaKey.orientation2color(R.parentBefore.meanOrientation));
    hold on; plot(R.parentBefore.boundary,'lineWidth',0.25,'lineColor','white'); hold off;
    title('Beta fragments | before mergeSimilar/mergeInclusions');
    saveFigure(outputFolder,'B3_Beta_BEFORE_merge.png');

    figure('Name','B4 Final beta IPF-Z','Color','white');
    plot(R.parentGrains,betaKey.orientation2color(R.parentGrains.meanOrientation));
    hold on; plot(R.parentGrains.boundary,'lineWidth',0.2,'lineColor','white'); hold off;
    title('Reconstructed beta (job.parentGrains) | IPF-Z');
    saveFigure(outputFolder,'B4_ParentBeta_IPFz.png');

    figure('Name','B5 Beta boundaries on alpha','Color','white');
    plot(alphaData,alphaKey.orientation2color(alphaData.orientations));
    hold on; plot(R.parentGrains.boundary,'lineWidth',0.1,'lineColor','black'); hold off;
    title('Measured alpha | inferred beta boundaries');
    saveFigure(outputFolder,'B5_BetaBoundaries_on_Alpha.png');

    if R.nRec > 0
        figure('Name','B6 Beta per pixel IPF-Z','Color','white');
        plot(R.rec,betaKey.orientation2color(R.rec.orientations));
        title('Reconstructed beta per pixel (converted alpha pixels only)');
        saveFigure(outputFolder,'B6_Beta_IPFz_PIXELS.png');

        figure('Name','B7 BOR misfit per pixel','Color','white');
        plot(R.rec,R.fit_deg); mtexColorbar; setColorRange([0 5]);
        hold on; plot(R.parentGrains.boundary,'lineWidth',0.5,'lineColor','white'); hold off;
        title('Measured alpha vs reconstructed beta misfit [deg]');
        saveFigure(outputFolder,'B7_Beta_PixelFit.png');
    end
catch ME
    warning('Failure in the section-6 figures: %s',ME.message);
end

%% 7. DIAGNOSTIC A: pixel correspondence check (ids and positions)
if R.nRec > 0
    alphaMatched = alphaData(ismember(alphaData.id,R.rec.id));
    [~,ia] = sort(alphaMatched.id);  [~,ib] = sort(R.rec.id);
    sameCount = length(alphaMatched) == R.nRec;
    sameId = sameCount && isequal(alphaMatched.id(ia),R.rec.id(ib));
    samePos = sameCount && max(abs(alphaMatched.x(ia)-R.rec.x(ib))) == 0 && ...
                           max(abs(alphaMatched.y(ia)-R.rec.y(ib))) == 0;
    fprintf('\n--- CORRESPONDENCE CHECK ---\n');
    fprintf('Same count: %d | same ids: %d | same x,y positions: %d\n', ...
        sameCount,sameId,samePos);
    if ~(sameCount && sameId && samePos)
        error('Alpha/beta correspondence inconsistent: do NOT use the PFs or the coverage.');
    end
else
    alphaMatched = alphaData([]);
end

%% 8. DIAGNOSTIC B: coverage (category map)
fprintf('\n--- RECONSTRUCTION COVERAGE ---\n');
fprintf('Indexed alpha (whole map):            %d px\n',R.nAlphaRaw);
fprintf('Alpha removed by minPixels:           %d px\n',R.nAlphaRaw-R.nAlphaUsed);
fprintf('Alpha with reconstructed beta:        %d px\n',R.nRec);
fprintf('Coverage (of indexed alpha):          %.2f %%\n',R.coverage_pct);
fprintf('Coverage (of alpha after cleanup):    %.2f %%\n',100*R.nRec/max(1,R.nAlphaUsed));

try
    ids = ebsdRaw.id;
    cat = ones(length(ebsdRaw),1);               % 1 not indexed
    cat(ismember(ids,alphaData.id)) = 2;         % 2 alpha NOT reconstructed
    cat(ismember(ids,R.rec.id))     = 3;         % 3 alpha reconstructed
    cat(ismember(ids,betaData.id))  = 4;         % 4 measured beta
    cmap = [0.85 0.85 0.85; 0.85 0.40 0.10; 0.15 0.45 0.80; 0.10 0.65 0.30];
    % Keep Coverage (13) in its OWN undocked window.  Do not let it
    % occupy / hide the later six-panel Pole Figures (12) in figure tabs.
    fig13Handle = figure('Name','C1 | Reconstruction coverage map', ...
        'NumberTitle','off','Color','white','WindowStyle','normal', ...
        'Visible','on','Tag','Ti64_Coverage_Fig13');
    figure(fig13Handle);
    plot(ebsdRaw,cat); colormap(gca,cmap); setColorRange([0.5 4.5]);
    hold on;
    hl = gobjects(4,1);
    for k = 1:4, hl(k) = patch(NaN,NaN,cmap(k,:)); end
    legend(hl,{'Not indexed','Alpha NOT reconstructed','Alpha reconstructed','Measured beta'}, ...
        'Location','southoutside','Orientation','horizontal');
    hold off;
    title(sprintf('Coverage: %.1f%% of indexed alpha',R.coverage_pct));
    saveFigure(outputFolder,'C1_Coverage_Map.png');
catch ME
    warning('Coverage map not generated: %s',ME.message);
end

fprintf('\nCoverage by size of the alpha child grain (pixel-weighted):\n');
disp(R.covByClass);
writetable(R.covByClass,fullfile(outputFolder,'Coverage_by_AlphaGrainSize.csv'));
fprintf('Median alpha ECD (area-weighted) = %.2f um -> %.1f px per diameter\n', ...
    R.alphaECD_areaW,R.alphaPxPerDiameter);
Tc = R.covByClass;
grp = {1:3,'< 3 um'; 4:5,'3-8 um'; 6:7,'>= 8 um'};
fprintf('Pooled coverage by alpha grain size:\n');
for g = 1:size(grp,1)
    ix = grp{g,1};
    fprintf('  alpha grains %-7s: %6d px (%4.1f%% of alpha) | reconstructed %.1f%%\n', ...
        grp{g,2},sum(Tc.AlphaPixels(ix)),100*sum(Tc.AlphaPixels(ix))/sum(Tc.AlphaPixels), ...
        100*sum(Tc.ReconstructedPixels(ix))/max(1,sum(Tc.AlphaPixels(ix))));
end
fprintf('(Low coverage in the coarse classes is expected if those grains are primary alpha.)\n');

%% 9. DIAGNOSTIC C: BOR quality
fprintf('\n--- BOR QUALITY (reconstructed pixels) ---\n');
if R.nRec > 0 && ~isnan(R.p50)
    fprintf('Median: %.2f deg | P90: %.2f deg | P99: %.2f deg\n',R.p50,R.p90,R.p99);
    fprintf('Misfit <= 2 deg: %.1f %% | <= 5 deg: %.1f %%\n',R.pct_le2,R.pct_le5);
    fprintf('Read TOGETHER with coverage (%.1f%%): low misfit with low coverage only\n',R.coverage_pct);
    fprintf('describes the easiest-to-reconstruct fraction.\n');
end

%% 10. DIAGNOSTIC D: fragmentation and beta grain size
% Beta grains are computed ONLY from alpha->beta pixels on the full map, so
% holes and edges count as boundary and isBoundary is meaningful.
S = R.frag;
fprintf('\n--- FRAGMENTATION (contiguous beta grains from reconstructed pixels) ---\n');
if S.n > 0
    fprintf('Total reconstructed beta grains: %d  (= %d independent orientations)\n',S.n,S.n);
    fprintf('Pixels per grain (measured): mean %.1f | median %.1f\n',S.meanPx,S.medPx);
    fprintf('1 px: %.1f%% | <= 3 px: %.1f%% | <= 10 px: %.1f%%\n',S.pct1,S.pctLE3,S.pctLE10);
    fprintf('Touching edge/holes: %d (%.0f%%) | interior: %d\n', ...
        nnz(S.edge),100*mean(S.edge),S.nInterior);
    fprintf('ECD [um]: mean(all) %.3f | mean(interior) %.3f | median(interior) %.3f | P90(all) %.3f\n', ...
        S.meanECD,S.meanECDint,S.medECDint,quantile(S.ecd,0.9));
    fprintf('Median ECD / step = %.1f  (1 px = %.2f steps; ~1-3 => pixel-dominated size)\n', ...
        median(S.ecd)/stepUm,2/sqrt(pi));
    fprintf('Fill ratio (measured px / area px) = %.2f  (<1: interior holes in the grains)\n',S.fill);
    fprintf('Interior-only statistics are biased LOW: large grains touch the map edge.\n');
    try
        figure('Name','C2 Fragmentation','Color','white');
        edgesPx = [0.5 1.5 2.5 3.5 5.5 10.5 20.5 50.5 100.5 500.5 Inf];
        cnt = histcounts(S.nPx,edgesPx);
        bar(cnt); set(gca,'XTickLabel',{'1','2','3','4-5','6-10','11-20','21-50','51-100','101-500','>500'});
        xlabel('Pixels per beta grain'); ylabel('Number of grains');
        title('Pixels per reconstructed beta grain');
        saveFigure(outputFolder,'C2_Fragmentation.png');

        figure('Name','C3 ECD histograms','Color','white');
        subplot(2,2,1); plotHist(S.ecd,'All');
        subplot(2,2,2); plotHist(S.ecd(~S.edge),'Interior');
        subplot(2,2,3); plotHist(S.ecd(S.edge),'Touching edge/holes');
        subplot(2,2,4); plotHist(S.ecd(S.nPx<=3),'Fragments <= 3 px');
        saveFigure(outputFolder,'C3_Beta_ECD_Histograms.png');
    catch ME
        warning('Histograms not generated: %s',ME.message);
    end
    Tg = table((1:S.n)',S.nPx(:),S.area(:),S.ecd(:),S.edge(:),S.nPx(:)<=3, ...
        'VariableNames',{'GrainIdx','Pixels','Area_um2','ECD_um','TouchesEdgeOrHole','Fragment_le3px'});
    writetable(Tg,fullfile(outputFolder,'ParentBeta_GrainMeasurements.csv'));
else
    warning('No contiguous beta grains for statistics.');
end
fprintf('\n(job.parentGrains: %d grains; mean ECD %.3f um; may include measured beta and\n', ...
    R.jobN,R.jobMeanECD);
fprintf(' is not comparable 1:1 with the grains above.)\n');

%% 11. Reliability
V = assessReliability(R,Q,mapMin);

%% 12. Pole figures: SIX-PANEL comparison + peak MRD and directions (v4.6)
% Critical fixes:
% 1. Individual PF failures (10 or 11) cannot prevent Figure D4.
% 2. Each of the SIX axes has an explicit MATLAB subplot parent. No nextAxis.
% 3. Figure D4 is TWO ROWS x THREE COLUMNS; alpha above, beta below.
% 4. One MRD scale for all three alpha PFs, a separate one for all beta PFs.
% 5. Figure D4 remains open and is saved to PNG AND editable FIG.
fig12Handle = [];
fig25Handle = [];
fig26Handle = [];
fig27Handle = [];
if R.nRec > 0
    hB = [Miller(1,1,0,CSbeta), Miller(2,0,0,CSbeta), Miller(1,1,1,CSbeta)];
    if strcmp(pfAlphaPopulation,'whole'), alphaPF = alphaData; else, alphaPF = alphaMatched; end

    % Figure D6. Failure here must NOT prevent Figure D4.
    try
        odfBetaPixel = calcDensity(R.rec.orientations,'halfwidth',pfHalfwidth,'silent');
        figure('Name','D6 Beta PF pixel-weighted','Color','white');
        plotPDF(odfBetaPixel,hB,'antipodal','contourf');
        mtexColorbar('title','mrd');
        saveFigure(outputFolder,'D6_ParentBeta_PF_PixelWeighted.png');
    catch ME10
        warning('Figure D6 could not be generated: %s',ME10.message);
    end

    % Figure D7. Independent from Figure D4.
    if S.n > 0
        try
            odfBetaArea = calcDensity(S.ori,'weights',S.area, ...
                'halfwidth',pfHalfwidth,'silent');
            figure('Name','D7 Beta PF grain-area-weighted','Color','white');
            plotPDF(odfBetaArea,hB,'antipodal','contourf');
            mtexColorbar('title','mrd');
            saveFigure(outputFolder,'D7_ParentBeta_PF_GrainAreaWeighted.png');
        catch ME11
            warning('Figure D7 could not be generated: %s',ME11.message);
        end
    end

    % Figure D4. Three measured-alpha PFs + three reconstructed-beta PFs.
    % SAME matched EBSD pixel IDs for the two orientation distributions.
    try
        fig12Handle = plotSixMatchedPFs( ...
            alphaPF.orientations,R.rec.orientations,CSalpha,CSbeta, ...
            pfHalfwidth,R.coverage_pct,S.n,outputFolder);
    catch ME12
        warning('Figure D4 SIX-PANEL plot failed: %s',ME12.message);
        fprintf(2,'FIGURE D4 ERROR DETAILS:\n%s\n', ...
            getReport(ME12,'extended','hyperlinks','off'));
    end

    % Figure D5. Add the ORIGINAL measured beta texture as an independent
    % three-row, NINE-PF comparison. Keep Figure D4 (paired 6 PFs) unchanged.
    % Row 2 uses ALL indexed beta points from the original CTF, not R.rec.
    % Do not compare pixel counts as if they were independent beta grains.
    try
        fig25Handle = plotNinePF_MeasuredAndReconstructedBeta( ...
            alphaPF.orientations,betaData.orientations, ...
            R.rec.orientations,CSalpha,CSbeta,pfHalfwidth, ...
            R.coverage_pct,outputFolder);
    catch ME25
        warning('Figure D5 NINE-PF plot failed: %s',ME25.message);
        fprintf(2,'FIGURE D5 ERROR DETAILS:\n%s\n', ...
            getReport(ME25,'extended','hyperlinks','off'));
    end

    if R.coverage_pct < 50
        warning(['The beta PF describes only %.1f%% of indexed alpha: NOT ', ...
                 'the beta texture of the whole sample.'],R.coverage_pct);
    end
    if S.n < 100
        warning(['The beta PF uses only %d reconstructed beta fragments; ', ...
                 'do not treat matched pixels as independent grain orientations.'],S.n);
    end
end

%% 13. Whole-map alpha PF check against Channel 5 (fig D1 + halfwidth scan)
% Fig. 28 is the whole-map measured alpha (all indexed alpha px, no matching):
% the same population Channel 5 uses for the alpha phase.  Figs 29/30 repeat fig 12.
try
    [fig28Handle,fig29Handle,fig30Handle] = plotThreeAlphaBetaPopulations( ...
        alphaData,alphaMatched,R.rec,CSalpha,CSbeta, ...
        pfHalfwidth,R.coverage_pct,showRedundantPF,outputFolder);
catch ME3pop
    warning('Whole-map alpha PF failed: %s',ME3pop.message);
    fprintf(2,'WHOLE-MAP ALPHA PF ERROR DETAILS:\n%s\n', ...
        getReport(ME3pop,'extended','hyperlinks','off'));
end
if runHalfwidthCheck
    try
        pfHalfwidthSweep(alphaData.orientations,CSalpha,halfwidthList_deg, ...
            ch5Ref0002,outputFolder);
    catch MEsw
        warning('Halfwidth check failed: %s',MEsw.message);
    end
end

% Figure D2/27: display the ORIGINAL measured beta phase alone.
% These figures are INDEPENDENT of the 9-PF comparison (Figure D5) and
% use only CTF pixels indexed as beta, never reconstructed pixels.
try
    [fig26Handle,fig27Handle] = plotMeasuredBetaPFandIPF( ...
        betaData,pfHalfwidth,showRedundantPF,outputFolder);
catch ME26
    warning('Measured-beta standalone PF/IPF figures failed: %s',ME26.message);
    fprintf(2,'MEASURED BETA FIGURES ERROR DETAILS:\n%s\n', ...
        getReport(ME26,'extended','hyperlinks','off'));
end

%% 13b. Scatter inverse pole figures (D8)
if runScatterPF && R.nRec > 0
    try
        plotScatterIPF(alphaData.orientations,CSalpha,betaData.orientations,CSbeta, ...
            R.rec.orientations,pfHalfwidth,scatterMaxPoints,outputFolder);
    catch MEsc
        warning('Scatter IPF figure failed: %s',MEsc.message);
    end
end

%% 14. PAIRED BOR VALIDATION: ALPHA (0001) vs NEAREST BETA {110}
% This analysis must not mix independent alpha and beta pixel populations.
% The pixel ID and XY identity test in section 7 is a required precondition.
% IMPORTANT: the nearest beta {110} pole is a MEMBER of the full family,
% not the full reconstructed {110} pole figure. Displaying only this member
% visualizes the plane-normal correspondence assumed by the BOR.
% parentEBSD.fit provides the *full* fitted Burgers-variant misorientation,
% whereas the plane-normal deviation cannot test the in-plane direction.
if runBORPairwiseCheck && R.nRec > 0
    try
        BOR = analyzePairedBOR(alphaMatched,R.rec,CSalpha,CSbeta, ...
            R.coverage_pct,BORopt,outputFolder);
    catch MEBOR
        warning('Paired BOR check failed: %s',MEBOR.message);
        fprintf('BOR details: %s\n',getReport(MEBOR,'extended','hyperlinks','off'));
    end
end

%% 15. SCHMID FACTOR: external load pointing LEFT (-X on EBSD map)
% Calculates SF for ALL measured alpha pixels, independently of BOR coverage.
% Reconstructed beta is a HISTORICAL parent, not today's alpha slip carrier.
% The Schmid factor predicts geometric resolved shear, NOT observed slip or
% a mechanical response. For genuine centrifugal loading, the full local
% stress tensor depends on specimen geometry, rotation and boundary conditions.
if runSchmidAnalysis
    try
        SFopt.ctfFile = fname;   % needed to link the macrozone (MTR) results
        SFresults = schmidAlphaAnalysis(alphaData,CSalpha,R,SFopt,outputFolder);
    catch ME
        warning('Schmid-factor analysis failed: %s',ME.message);
        fprintf('Schmid detail: %s\n',getReport(ME,'extended','hyperlinks','off'));
    end
end

%% 16. Comparison with your original parameters
cmpRows = metricsRow(R,P,stepUm,'main');
if runBaselineComparison
    Rb = runReconstruction(ebsdRaw,alphaPhase,betaPhase,Pbase,stepUm);
    cmpRows = [cmpRows; metricsRow(Rb,Pbase,stepUm,'baseline_original')];
    fprintf('\n--- COMPARISON: new vs original parameters ---\n');
    disp(cmpRows);
    writetable(cmpRows,fullfile(outputFolder,'Comparison_main_vs_baseline.csv'));
end

%% 17. Parameter sensitivity (optional)
if runSensitivity
    list = buildSensitivityList(P,sensMode);
    Tsens = table();
    for i = 1:numel(list)
        fprintf('Sensitivity %d/%d ...\n',i,numel(list));
        try
            Ri = runReconstruction(ebsdRaw,alphaPhase,betaPhase,list{i},stepUm);
        catch ME
            Ri = struct('ok',false,'reason',ME.message);
        end
        Tsens = [Tsens; metricsRow(Ri,list{i},stepUm,sprintf('run%03d',i))]; %#ok<AGROW>
        writetable(Tsens,fullfile(outputFolder,'Sensitivity_parameters.csv'));
    end
    disp(Tsens);
end

%% 18. Convergence versus step size (optional)
% Subsamples the SAME map and repeats the reconstruction. If the beta size is
% physical, the median ECD [um] should plateau as the step grows; if it grows
% in proportion to the step, it is a pixel artefact.
if runConvergence
    Tconv = table();
    for k = convFactors
        fprintf('Convergence: factor %d ...\n',k);
        try
            if k == 1, ek = ebsdRaw; else, ek = downsampleEBSD(ebsdRaw,k); end
            sk = estimateStep(ek);
            Rk = runReconstruction(ek,alphaPhase,betaPhase,P,sk);
        catch ME
            Rk = struct('ok',false,'reason',ME.message); sk = stepUm*k;
        end
        Tconv = [Tconv; metricsRow(Rk,P,sk,sprintf('k%d',k))]; %#ok<AGROW>
    end
    disp(Tconv);
    writetable(Tconv,fullfile(outputFolder,'Convergence_vs_StepSize.csv'));
    try
        figure('Name','G1 Convergence','Color','white');
        loglog(Tconv.Step_um,Tconv.MedianECD_um,'o-',Tconv.Step_um,Tconv.MeanECD_um,'s-'); hold on;
        loglog(Tconv.Step_um,Tconv.Step_um*Tconv.MedianECD_um(1)/Tconv.Step_um(1),'k--');
        hold off; grid on;
        xlabel('Step size [um]'); ylabel('Beta ECD [um]');
        legend('Median','Mean','Slope 1 (size = k x pixel)','Location','northwest');
        title('If it follows the dashed line, the size is a pixel artefact');
        saveFigure(outputFolder,'G1_Convergence.png');
    catch ME
        warning('Convergence figure not generated: %s',ME.message);
    end
end

%% 19. Save
settings = struct('P',P,'Pbase',Pbase,'Q',Q,'inputFile',fname,'stepUm',stepUm, ...
    'roiRect',roiRect,'downsampleFactor',downsampleFactor,'reliability',V, ...
    'runSchmidAnalysis',runSchmidAnalysis,'SFopt',SFopt, ...
    'runBORPairwiseCheck',runBORPairwiseCheck,'BORopt',BORopt);
Rsave = rmfield(R,intersect(fieldnames(R),{'job'}));
save(fullfile(outputFolder,'ParentBeta_v4_BOR_Schmid_Workspace.mat'),'Rsave','settings','-v7.3');
fprintf('\nDONE. Results in:\n%s\n',outputFolder);

%% 20. Display Figure D4 separately (never underneath Coverage Figure 13)
% MATLAB figure *numbers* do not necessarily match our saved filenames:
% rely on the exact figure HANDLE / Tag rather than gcf or figure(12).
if exist('fig13Handle','var') && isgraphics(fig13Handle,'figure')
    set(fig13Handle,'WindowStyle','normal','Visible','on');
end
if ~isempty(fig12Handle) && isgraphics(fig12Handle,'figure')
    set(fig12Handle,'WindowStyle','normal','Visible','on');
    figure(fig12Handle); drawnow;
    try, shg; catch, end
    fprintf('\nFIGURE D4 OPEN IN ITS OWN WINDOW (independent of Coverage 13).\n');
    fprintf('MATLAB figure number = %g; name = %s\n', ...
        get(fig12Handle,'Number'),get(fig12Handle,'Name'));
    fprintf('Figure D4 axes found: %d (expected 6).\n', ...
        numel(findall(fig12Handle,'Type','axes')));
    fprintf('Saved in: %s\n',outputFolder);
elseif exist(fullfile(outputFolder,'D4_PF_Alpha_vs_Beta.fig'),'file') == 2
    warning('Figure D4 window closed; reopening saved six-panel FIG.');
    fig12Handle = openfig(fullfile(outputFolder, ...
        'D4_PF_Alpha_vs_Beta.fig'),'visible');
    set(fig12Handle,'WindowStyle','normal','Tag','Ti64_SixPF_Fig12');
    figure(fig12Handle); drawnow;
else
    warning('Figure D4 could not be generated; inspect FIGURE D4 ERROR DETAILS above.');
end

% Figure D5 is a different MATLAB window, not an overlay on Figures 12/13.
if ~isempty(fig25Handle) && isgraphics(fig25Handle,'figure')
    set(fig25Handle,'WindowStyle','normal','Visible','on');
    fprintf('Figure D5 (NINE PF) is open independently from Figures 12 and 13.\n');
    figure(fig25Handle); drawnow;
end
% Measured-beta IPFs are brought to the foreground LAST: each figure is
% in its own independent MATLAB window (no Figure 13 overlap).
if ~isempty(fig26Handle) && isgraphics(fig26Handle,'figure')
    set(fig26Handle,'WindowStyle','normal','Visible','on');
end
if ~isempty(fig27Handle) && isgraphics(fig27Handle,'figure')
    set(fig27Handle,'WindowStyle','normal','Visible','on');
    figure(fig27Handle); drawnow;
end

%% ======================= LOCAL FUNCTIONS =================================

function R = runReconstruction(ebsdRaw,alphaPhase,betaPhase,P,step)
% Full reconstruction + metrics. Produces no figures.
R = struct('ok',false,'reason','');
alphaRaw = ebsdRaw(alphaPhase);
betaRaw  = ebsdRaw(betaPhase);
R.step = step;  R.nAlphaRaw = length(alphaRaw);
alphaId = alphaRaw.phaseId(1);

ebsd = ebsdRaw('indexed');
[grains,ebsd.grainId] = calcGrains(ebsd,'threshold',P.childThreshold);
if P.minPixels > 1
    tiny = grains(grains.grainSize < P.minPixels & grains.phaseId == alphaId);
    ebsd(tiny) = [];
    [grains,ebsd.grainId] = calcGrains(ebsd,'threshold',P.childThreshold);
end
R.nAlphaUsed = nnz(ismember(ebsd.id,alphaRaw.id));

p2c0 = orientation.Burgers(betaRaw.CS,alphaRaw.CS);
job = parentGrainReconstructor(ebsd,grains,p2c0);
job.p2c = p2c0;
R.gbFit0 = finiteVec(job.calcGBFit./degree);
if P.refineBOR && numel(R.gbFit0) >= 20
    job.calcParent2Child('quantile',P.refQuantile);
    if angle(p2c0,job.p2c)./degree > 5, job.p2c = p2c0; end
end
R.gbFit = finiteVec(job.calcGBFit./degree);

args = {'threshold',P.graphThreshold};
if ~isempty(P.graphTolerance), args = [args,{'tolerance',P.graphTolerance}]; end
job.calcVariantGraph(args{:});
job.clusterVariantGraph('numIter',P.graphIterations);
job.calcParentFromVote('minProb',P.minProb);

nChild = nnz(job.grainsPrior.phaseId == job.childPhaseId);
R.alphaGrainsAssigned_pct = 100*nnz(job.isTransformed)/max(1,nChild);
if nnz(job.isTransformed) == 0
    R.reason = 'no alpha grain assigned after voting'; return
end

R.parentBefore = job.parentGrains;
job.mergeSimilar('threshold',P.mergeThreshold);
job.mergeInclusions('maxSize',P.maxInclusionPixels);
R.parentGrains = job.parentGrains;
if isempty(R.parentGrains), R.reason = 'parentGrains is empty'; return, end
R.job = job;
parentEBSD = job.calcParentEBSD;

% Beta pixels that come from alpha (excludes measured beta)
recAll = parentEBSD(betaPhase);
if isempty(recAll)
    R.rec = recAll;
else
    R.rec = recAll(ismember(recAll.id,alphaRaw.id));
end
R.nRec = length(R.rec);
R.coverage_pct = 100*R.nRec/R.nAlphaRaw;

% BOR quality
R.fit_deg = []; R.p50 = NaN; R.p90 = NaN; R.p99 = NaN; R.pct_le2 = NaN; R.pct_le5 = NaN;
if R.nRec > 0
    f = R.rec.fit./degree;  R.fit_deg = f;
    ff = f(isfinite(f));
    if ~isempty(ff)
        q = quantile(ff,[0.5 0.9 0.99]);
        R.p50 = q(1); R.p90 = q(2); R.p99 = q(3);
        R.pct_le2 = 100*mean(ff <= 2);  R.pct_le5 = 100*mean(ff <= 5);
    end
end

% Alpha child statistics and coverage by alpha grain size
ecd = 2*sqrt(area(grains)./pi);
[~,loc] = ismember(ebsd.grainId,grains.id);
pixECD = nan(length(ebsd),1);  pixECD(loc>0) = ecd(loc(loc>0));
isA = ismember(ebsd.id,alphaRaw.id);
isR = false(length(ebsd),1);
if R.nRec > 0, isR = ismember(ebsd.id,R.rec.id); end
R.alphaECD_areaW = median(pixECD(isA & isfinite(pixECD)));
R.alphaPxPerDiameter = R.alphaECD_areaW/step;
edgesE = [0 1 2 3 5 8 12 Inf];
nCls = numel(edgesE)-1;
nPxCls = zeros(nCls,1); nRecCls = zeros(nCls,1);
for c = 1:nCls
    m = isA & pixECD >= edgesE(c) & pixECD < edgesE(c+1);
    nPxCls(c) = nnz(m);  nRecCls(c) = nnz(m & isR);
end
lbl = arrayfun(@(a,b) sprintf('%g-%g um',a,b),edgesE(1:end-1),edgesE(2:end),'UniformOutput',false)';
R.covByClass = table(lbl,nPxCls,nRecCls,100*nRecCls./max(1,nPxCls), ...
    'VariableNames',{'AlphaGrainECD','AlphaPixels','ReconstructedPixels','Coverage_pct'});

% Beta fragments (from reconstructed pixels, on the full map)
R.frag = fragmentStats(parentEBSD,betaRaw,R.rec,P.mergeThreshold,step);

% job.parentGrains (reference only)
R.jobN = length(R.parentGrains);
R.jobMeanECD = mean(2*sqrt(area(R.parentGrains)./pi));
R.ok = true;
end

function S = fragmentStats(parentEBSD,betaRaw,rec,thr,step)
S = struct('n',0,'meanPx',NaN,'medPx',NaN,'pct1',NaN,'pctLE3',NaN,'pctLE10',NaN, ...
    'meanECD',NaN,'meanECDint',NaN,'medECDint',NaN,'nInterior',0,'fill',NaN, ...
    'nPx',[],'ecd',[],'area',[],'edge',logical([]),'ori',[]);
if isempty(rec), return, end
pe = parentEBSD('indexed');
pe(ismember(pe.id,betaRaw.id)) = [];      % remove MEASURED beta
try
    bg = calcGrains(pe,'threshold',thr);
    bg = bg(bg.phaseId == rec.phaseId(1));
    A = area(bg);  A = A(:);
    nPx = bg.grainSize;  nPx = nPx(:);    % measured pixels (excludes interior holes)
    edge = logical(isBoundary(bg));  edge = edge(:);
    ori = bg.meanOrientation;
catch
    return
end
ecd = 2*sqrt(A./pi);
S.n = numel(A);  S.nPx = nPx;  S.area = A;  S.ecd = ecd;  S.edge = edge;  S.ori = ori;
S.fill = sum(nPx)/sum(A./step^2);
S.meanPx = mean(nPx);  S.medPx = median(nPx);
S.pct1 = 100*mean(nPx == 1);  S.pctLE3 = 100*mean(nPx <= 3);  S.pctLE10 = 100*mean(nPx <= 10);
S.meanECD = mean(ecd);
S.nInterior = nnz(~edge);
if S.nInterior > 0
    S.meanECDint = mean(ecd(~edge));  S.medECDint = median(ecd(~edge));
end
end

function V = assessReliability(R,Q,mapMin)
% Heuristic verdict (criteria in Q, adjustable). Prints PASS/FAIL.
V = struct('pass',struct());
S = R.frag;
if S.n == 0
    V.reliableGrainSize = false;
    fprintf('\n=== RELIABILITY VERDICT: no beta grains to assess ===\n');
    return
end
ecdP90 = quantile(S.ecd,0.9);
pctEdge = 100*mean(S.edge);
pctImpl = 100*mean(S.ecd < Q.minPlausibleBetaECD_um);
V.pass.sampling   = R.alphaPxPerDiameter >= Q.minPxPerAlphaDiameter;
V.pass.bor        = isfinite(R.p90) && R.p50 <= Q.maxMedianBOR_deg && R.p90 <= Q.maxP90BOR_deg;
V.pass.fragments  = S.medPx >= Q.minMedianFragmentPx && S.pctLE3 <= Q.maxPctFragLE3;
V.pass.plausible  = pctImpl <= Q.maxPctImplausible;
V.pass.statistics = S.nInterior >= Q.minInteriorGrains;
V.pass.fov        = mapMin >= Q.minFovOverBetaECD*ecdP90;
V.pass.edge       = pctEdge <= Q.maxPctEdge;
V.pass.coverage   = R.coverage_pct >= 50;
fprintf('\n=== RELIABILITY VERDICT (heuristic) ===\n');
pf = @(b) ternary(b,'PASS','FAIL');
fprintf('Sampling (%.1f px per alpha diameter, min %g):         %s\n',R.alphaPxPerDiameter,Q.minPxPerAlphaDiameter,pf(V.pass.sampling));
fprintf('BOR (P50 %.2f, P90 %.2f deg):                          %s\n',R.p50,R.p90,pf(V.pass.bor));
fprintf('Fragments (median %.0f px, %.0f%% <= 3 px):            %s\n',S.medPx,S.pctLE3,pf(V.pass.fragments));
fprintf('Plausible size (%.0f%% of grains < %g um):             %s\n',pctImpl,Q.minPlausibleBetaECD_um,pf(V.pass.plausible));
fprintf('Statistics (%d interior grains, min %d):               %s\n',S.nInterior,Q.minInteriorGrains,pf(V.pass.statistics));
fprintf('Field of view (side %.0f um vs %g x P90 ECD = %.0f um): %s\n',mapMin,Q.minFovOverBetaECD,Q.minFovOverBetaECD*ecdP90,pf(V.pass.fov));
fprintf('Truncation (%.0f%% of grains touch edge/holes, max %g%%): %s\n',pctEdge,Q.maxPctEdge,pf(V.pass.edge));
fprintf('Coverage >= 50%% (%.1f%%) [informative]:                %s\n',R.coverage_pct,pf(V.pass.coverage));
core = V.pass.sampling && V.pass.bor && V.pass.fragments && V.pass.plausible && ...
       V.pass.statistics && V.pass.fov && V.pass.edge;
V.reliableGrainSize = core;
if core
    fprintf('>> Beta grain size: BASIC CRITERIA MET (still validate morphology and primary alpha).\n');
else
    fprintf('>> NOT RELIABLE as prior-beta grain size: these are BOR-compatible regions,\n');
    fprintf('   not identified prior-beta grains. Do not use these ECDs in conclusions.\n');
    fprintf('   (The local orientations can still be informative; the SIZE is not.)\n');
end
if ~V.pass.coverage
    fprintf('>> Reconstruction covers a minority of the alpha: beta PFs are NOT global texture.\n');
end
end

function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end

function plotHist(x,ttl)
if isempty(x)
    text(0.5,0.5,'none','HorizontalAlignment','center'); axis off; title(ttl);
else
    histogram(x,'BinWidth',1); xlabel('ECD [um]'); ylabel('Grains'); title(ttl);
end
end

function row = metricsRow(R,P,step,name)
nan1 = NaN;
if isempty(P.graphTolerance), tol = NaN; else, tol = P.graphTolerance/degree; end
base = {name,step,P.childThreshold/degree,P.graphThreshold/degree,tol,P.minProb, ...
        P.minPixels,P.mergeThreshold/degree};
if isfield(R,'ok') && R.ok && R.frag.n > 0
    S = R.frag;
    v = {S.n,R.coverage_pct,R.p50,R.p90,S.meanECD,median(S.ecd),median(S.ecd)/step, ...
         S.pct1,S.pctLE3,S.pctLE10,R.alphaPxPerDiameter,100*mean(S.edge),S.nInterior};
else
    v = {0,nan1,nan1,nan1,nan1,nan1,nan1,nan1,nan1,nan1,nan1,nan1,nan1};
    if isfield(R,'ok') && R.ok, v{2} = R.coverage_pct; end
end
row = cell2table([base v],'VariableNames',{'Run','Step_um','ChildThreshold_deg', ...
    'BORGraphThreshold_deg','BORGraphTolerance_deg','MinVoteProb','MinPixels', ...
    'ParentMerge_deg','NBetaGrains','Coverage_pct','MedianBOR_deg','P90BOR_deg', ...
    'MeanECD_um','MedianECD_um','MedianECD_over_step','Pct_1px','Pct_le3px', ...
    'Pct_le10px','AlphaPxPerDiameter','Pct_TouchEdge','N_interior'});
end

function list = buildSensitivityList(P0,mode)
child = [1.5 2 3]*degree;  bor = [1.5 2 3]*degree;  mp = [0.30 0.50];
mpx = [1 3 5];  mrg = [3 5]*degree;
list = {};
add = @(c,b,m,x,g) setfields(P0,c,b,m,x,g);
if strcmp(mode,'fullGrid')
    for c = child, for b = bor, for m = mp, for x = mpx, for g = mrg
        list{end+1} = add(c,b,m,x,g); %#ok<AGROW>
    end, end, end, end, end
else
    list{end+1} = P0;
    for c = child, list{end+1} = add(c,P0.graphThreshold,P0.minProb,P0.minPixels,P0.mergeThreshold); end %#ok<AGROW>
    for b = bor,   list{end+1} = add(P0.childThreshold,b,P0.minProb,P0.minPixels,P0.mergeThreshold); end %#ok<AGROW>
    for m = mp,    list{end+1} = add(P0.childThreshold,P0.graphThreshold,m,P0.minPixels,P0.mergeThreshold); end %#ok<AGROW>
    for x = mpx,   list{end+1} = add(P0.childThreshold,P0.graphThreshold,P0.minProb,x,P0.mergeThreshold); end %#ok<AGROW>
    for g = mrg,   list{end+1} = add(P0.childThreshold,P0.graphThreshold,P0.minProb,P0.minPixels,g); end %#ok<AGROW>
end
end

function P = setfields(P,c,b,m,x,g)
P.childThreshold = c;  P.graphThreshold = b;  P.graphTolerance = b;  % tolerance = threshold
P.minProb = m;  P.minPixels = x;  P.mergeThreshold = g;
end

function s = estimateStep(ebsd)
ux = unique(round(ebsd.x,6));
d = diff(ux);  d = d(d > 0);
if isempty(d)
    uy = unique(round(ebsd.y,6)); d = diff(uy); d = d(d > 0);
end
s = min(d);
end

function e2 = downsampleEBSD(e,k)
% Regular-grid subsampling: every k-th pixel in x and y.
e2 = gridify(e);
e2 = e2(1:k:end,1:k:end);
end

function v = finiteVec(v)
v = v(:);  v = v(isfinite(v));
end

function d = ch5EquivDeg(hw)
% MTEX halfwidth [rad] -> equivalent Channel 5 "half width" [deg].
% Empirical factor 2 (SAMPLE015, Channel 5 at 10 and 20 deg). Not Oxford documentation.
d = 2*hw/degree;
end

function T = pfHalfwidthSweep(ori,cs,hwList_deg,ch5Ref,outFolder)
% Maximum of the alpha PFs vs MTEX halfwidth. Use it to pick the halfwidth that
% reproduces a Channel 5 maximum (ch5Ref = Channel 5 {0002} max, optional).
h = [Miller(0,0,0,2,cs),Miller(1,0,-1,0,cs),Miller(1,1,-2,0,cs)];
fam = {'{0002}','{10-10}','{11-20}'};
rows = cell(0,7);
for i = 1:numel(hwList_deg)
    hw = hwList_deg(i)*degree;
    odf = calcDensity(ori,'halfwidth',hw,'silent');
    for k = 1:3
        [pk,~,tilt,az] = pfMaxWithDirection(odf,h(k));
        rows(end+1,:) = {fam{k},hwList_deg(i),ch5EquivDeg(hw),pk,tilt,az,numel(ori)}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'PoleFamily','MTEX_halfwidth_deg', ...
    'Channel5_equiv_halfwidth_deg','MaxMRD','TiltFromZ_deg','AzimuthFromX_deg','EBSD_Pixels'});
writetable(T,fullfile(outFolder,'D0_PF_Halfwidth_Check.csv'));
fprintf('\n===== WHOLE-MAP ALPHA: PF MAXIMUM vs KERNEL =====\n');
disp(T(:,1:4));
if ~isempty(ch5Ref)
    k = strcmp(T.PoleFamily,'{0002}');
    [v,ix] = sort(T.MaxMRD(k));  hwv = T.MTEX_halfwidth_deg(k);  hwv = hwv(ix);
    for r = ch5Ref(:)'
        heq = exp(interp1(log(v),log(hwv),log(r),'linear'));
        if isnan(heq)
            warning('Channel 5 value %.2f is outside the scanned range.',r);
        else
            fprintf('Channel 5 {0002} = %.2f  <->  MTEX halfwidth ~%.2f deg (Channel 5 "%.1f deg")\n', ...
                r,heq,2*heq);
        end
    end
end
end

function figH = plotScatterIPF(oriA,csA,oriBm,csB,oriBr,halfwidth,maxPts,outFolder)
% D8: SCATTER inverse pole figure for the specimen Z direction (IPF-Z), one panel per
% phase: whole-map alpha / measured beta / reconstructed beta. One dot per orientation,
% in the fundamental sector of each phase. Large populations are randomly subsampled
% to maxPts dots (display only).
    oris = {oriA,oriBm,oriBr};
    names = {'Alpha (whole map)','Beta measured','Beta reconstructed'};
    figH = figure('Name','D8 | SCATTER IPF-Z - alpha / beta / reconstructed beta', ...
        'NumberTitle','off','Color','white','Visible','on','WindowStyle','normal', ...
        'Tag','Ti64_ScatterIPF_D8','Units','pixels','Position',[100 120 1700 560]);
    rng(11);
    for r = 1:3
        ori = oris{r};
        if isempty(ori), continue; end
        n = numel(ori);
        if n > maxPts, ori = ori(randperm(n,maxPts)); end
        ax = subplot(1,3,r,'Parent',figH);
        try
            plotIPDF(ori,zvector,'MarkerSize',3,'parent',ax);
        catch ME
            warning('Scatter IPF-Z failed (%s): %s',names{r},ME.message);
        end
        title(ax,sprintf('%s | IPF-Z | N = %d (%d dots)',names{r},n,numel(ori)), ...
            'FontSize',11,'FontWeight','bold','Interpreter','none');
    end
    drawnow;
    exportgraphics(figH,fullfile(outFolder,'D8_ScatterIPFZ_Alpha_Beta_BetaRec.png'),'Resolution',300);
    savefig(figH,fullfile(outFolder,'D8_ScatterIPFZ_Alpha_Beta_BetaRec.fig'));
end

function T = plotPhaseMap(ebsdAll,alphaData,betaData,sampleName,stepUm,outFolder)
% A2: PHASE MAP of the whole scan: alpha / beta / zero solution (not indexed).
% Zero solution = points that the indexing software could not solve (phase 0).
% Also prints and saves the area fraction of each class (of ALL scan points).
    ids = ebsdAll.id;
    cls = ones(length(ebsdAll),1);              % 1 zero solution
    cls(ismember(ids,alphaData.id)) = 2;        % 2 alpha
    cls(ismember(ids,betaData.id))  = 3;        % 3 beta (measured)
    cmap = [0.05 0.05 0.05; 0.15 0.45 0.80; 0.85 0.40 0.10];
    names = {'Zero solution','Alpha (Ti-Hex)','Beta (Titanium cubic)'};
    figH = figure('Name','A2 | Phase map: alpha / beta / zero solution', ...
        'NumberTitle','off','Color','white','Visible','on','WindowStyle','normal', ...
        'Tag','Ti64_PhaseMap_A2');
    plot(ebsdAll,cls); colormap(gca,cmap); setColorRange([0.5 3.5]);
    hold on;
    hl = gobjects(3,1);
    for k = 1:3, hl(k) = patch(NaN,NaN,cmap(k,:)); end
    legend(hl,names,'Location','southoutside','Orientation','horizontal');
    hold off;
    n = numel(cls);
    cnt = accumarray(cls,1,[3 1]);
    pct = 100*cnt/n;
    title(sprintf('%s | zero solution %.1f%% | alpha %.1f%% | beta %.2f%%', ...
        sampleName,pct(1),pct(2),pct(3)),'Interpreter','none');
    T = table(names(:),cnt,pct,cnt*stepUm^2, ...
        'VariableNames',{'Class','Points','Percent_of_scan','Area_um2'});
    disp(T);
    writetable(T,fullfile(outFolder,'A2_PhaseFractions.csv'));
    saveFigure(outFolder,'A2_PhaseMap.png');
    try, savefig(figH,fullfile(outFolder,'A2_PhaseMap.fig')); catch, end
end

function saveFigure(folder,name)
exportgraphics(gcf,fullfile(folder,name),'Resolution',300);
end



%% ---------- PAIRED BOR / POLE-FIGURE HELPERS ------------------------------
function B = analyzePairedBOR(alphaMatched,rec,CSalpha,CSbeta,coveragePct,opt,folder)
% The alpha / reconstructed-beta pair must have the SAME EBSD pixel id.
% The two crystal phases are never compared via unrelated unpaired PF maxima.
[idA,ia] = sort(alphaMatched.id(:));
[idB,ib] = sort(rec.id(:));
if numel(idA) ~= numel(idB) || ~isequal(idA,idB)
    error('BOR comparison aborted: measured alpha and reconstructed beta ids differ.');
end
A = alphaMatched(ia);
Brec = rec(ib);
if any(A.x(:) ~= Brec.x(:) | A.y(:) ~= Brec.y(:))
    error('BOR comparison aborted: matched ids do not have identical XY positions.');
end
n = numel(idA);

% Burgers: (0001)alpha || {110}beta.
% In hexagonal Miller notation this is the alpha c-axis [0001].
% Cubic [110] is normal to the (110) plane (this relies on cubic symmetry).
% The full OR ALSO constrains <11-20>alpha || <111>beta; see rec.fit.
alphaC = Miller(0,0,0,1,CSalpha,'UVTW');
beta110 = Miller(1,1,0,CSbeta,'uvw');
beta110Family = symmetrise(beta110,'unique','antipodal');
nb = length(beta110Family);
if nb < 3, error('Unexpected beta {110} family size: %d.',nb); end

% Each alpha orientation is compared with every symmetrically equivalent
% beta {110} axis, and the nearest one is identified. The antipodal flag
% means normals v and -v represent the same physical lattice plane.
minPlaneDeg = inf(n,1);
nearestPole = zeros(n,1,'uint8');
for startIdx=1:opt.chunkSize:n
    ii = startIdx:min(startIdx+opt.chunkSize-1,n);
    cSpec = A.orientations(ii) * alphaC;
    bo = Brec.orientations(ii);
    localMin = inf(numel(ii),1);
    localIdx = zeros(numel(ii),1,'uint8');
    for kk=1:nb
        bSpec = bo * beta110Family(kk);
        d = angle(cSpec,bSpec,'antipodal') ./ degree;
        d = d(:);
        take = d < localMin;
        localMin(take) = d(take);
        localIdx(take) = uint8(kk);
    end
    minPlaneDeg(ii) = localMin;
    nearestPole(ii) = localIdx;
end

% MTEX's reconstructed parentEBSD.fit evaluates the COMPLETE Burgers OR
% (the nearest compatible child variant), not only the basal-plane normal.
% It is a reconstruction-based metric, not independent experimental proof.
fullBORfit_deg = Brec.fit(:) ./ degree;
validPlane = isfinite(minPlaneDeg);
validFull = isfinite(fullBORfit_deg);
if ~any(validPlane) || ~any(validFull)
    error('No finite BOR metrics were obtained.');
end
qp = quantile(minPlaneDeg(validPlane),[0.5 0.9 0.99]);
qf = quantile(fullBORfit_deg(validFull),[0.5 0.9 0.99]);
tol = opt.toleranceDeg(:)';

B = struct();
B.id = idA;
B.x_um = A.x(:);
B.y_um = A.y(:);
B.PlaneNormalDeviation_deg = minPlaneDeg;
B.FullBORVariantFit_deg = fullBORfit_deg;
B.ClosestBeta110PoleIndex = nearestPole;
B.nPairs = n;
B.beta110FamilySize = nb;
B.coverage_pct = coveragePct;
B.planeQuantiles_deg = qp;
B.fullQuantiles_deg = qf;

fprintf('\n--- PAIRED BOR VALIDATION: alpha [0001] vs beta {110} ---\n');
fprintf('Matched alpha/beta pixels: %d (%.2f%% of indexed alpha)\n',n,coveragePct);
fprintf('Basal / beta-{110} normal deviation: median %.3f, P90 %.3f, P99 %.3f deg\n',qp);
fprintf('Full BOR variant fit (MTEX):          median %.3f, P90 %.3f, P99 %.3f deg\n',qf);
for k=1:numel(tol)
    fprintf('Fraction <= %g deg: plane %.2f%% | full BOR %.2f%%\n',tol(k), ...
        100*mean(minPlaneDeg(validPlane)<=tol(k)), ...
        100*mean(fullBORfit_deg(validFull)<=tol(k)));
end
pairedValid = validPlane & validFull;
inconsistent = pairedValid & (minPlaneDeg > fullBORfit_deg + 1.0);
fprintf('Plane normal > full BOR fit by >1 deg: %d of %d paired finite pixels\n', ...
    nnz(inconsistent),nnz(pairedValid));
if 100*nnz(inconsistent)/max(1,nnz(pairedValid)) > 1
    warning('Paired plane vs full-BOR mismatch: check EBSD frames / Miller notation.');
end
if coveragePct < 50
    warning('Paired BOR metrics describe only %.1f%% of indexed alpha, not the whole map.',coveragePct);
end
fprintf(['NOTE: low fit is EXPECTED when a BOR-constrained reconstruction succeeds.\n' ...
 'It is not an independent validation of the real parent beta grains.\n']);

% Machine-readable outputs, without duplicating a large pixel table by default.
Tsum = table(n,coveragePct,nb,qp(1),qp(2),qp(3),qf(1),qf(2),qf(3), ...
    100*mean(minPlaneDeg(validPlane)<=2),100*mean(minPlaneDeg(validPlane)<=5), ...
    100*mean(fullBORfit_deg(validFull)<=2),100*mean(fullBORfit_deg(validFull)<=5), ...
    100*nnz(inconsistent)/max(1,nnz(pairedValid)), ...
    'VariableNames',{'MatchedPixels','AlphaCoverage_pct','Beta110FamilySize', ...
    'PlaneMedian_deg','PlaneP90_deg','PlaneP99_deg', ...
    'FullBORMedian_deg','FullBORP90_deg','FullBORP99_deg', ...
    'PlaneLE2_pct','PlaneLE5_pct','FullBORLE2_pct','FullBORLE5_pct', ...
    'PlaneExceedsFullBy1deg_pct'});
B.summary = Tsum;
writetable(Tsum,fullfile(folder,'E0_BOR_Pairwise_Summary.csv'));
save(fullfile(folder,'E0_BOR_Pairwise_AllPixels.mat'),'B','-v7.3');
if opt.exportPixelCSV
    Tout = table(idA,A.x(:),A.y(:),minPlaneDeg,fullBORfit_deg,double(nearestPole), ...
        'VariableNames',{'EBSD_ID','X_um','Y_um','PlaneNormalDeviation_deg', ...
        'FullBORVariantFit_deg','SelectedBeta110FamilyIndex'});
    writetable(Tout,fullfile(folder,'E0_BOR_Pairwise_AllPixels.csv'));
end

% Fig. 18: scatter overlay of CORRESPONDING alpha c-axes and beta poles.
% We deliberately choose only the nearest beta {110} pole to each alpha
% basal normal: a full beta {110} PF would contain ALL equivalent poles.
% Both sets use a common specimen frame and equal-area upper-hemisphere
% projection; each axis is independently folded through its antipode.
try
    np = min(n,opt.maxOverlayPairs);
    oldRng = rng;
    rng(opt.randomSeed,'twister');
    pick = randperm(n,np).';
    rng(oldRng);
    cp = A.orientations(pick) * alphaC;
    [ax,ay] = borLambertUpper(cp);
    bx = nan(np,1); by = nan(np,1);
    for kk=1:nb
        sel = nearestPole(pick)==kk;
        if any(sel)
            bp = Brec.orientations(pick(sel)) * beta110Family(kk);
            [bx(sel),by(sel)] = borLambertUpper(bp);
        end
    end
    figure('Name','E1 Alpha basal / nearest beta110 PF overlay','Color','white');
    hold on;
    t = linspace(0,2*pi,512);
    plot(cos(t),sin(t),'k-','LineWidth',1.1,'HandleVisibility','off');
    plot([-1 1],[0 0],':','Color',[.80 .80 .80],'HandleVisibility','off');
    plot([0 0],[-1 1],':','Color',[.80 .80 .80],'HandleVisibility','off');
    nLines = min(opt.maxPairLines,np);
    if nLines > 0
        lineIdx = round(linspace(1,np,nLines));
        for j=lineIdx(:)'
            % Projection boundary antipodal folding can place a valid
            % pair on opposite sides of the disk. Do not connect such points.
            if hypot(ax(j)-bx(j),ay(j)-by(j)) < 0.25
                plot([ax(j),bx(j)],[ay(j),by(j)],'-', ...
                    'Color',[.72 .72 .72],'LineWidth',0.55,'HandleVisibility','off');
            end
        end
    end
    ha = plot(ax,ay,'.','Color',[.10 .33 .85],'MarkerSize',10);
    hb = plot(bx,by,'x','Color',[.92 .22 .12],'MarkerSize',5,'LineWidth',0.75);
    axis equal; axis([-1.08 1.08 -1.08 1.08]); axis off;
    text(0,1.055,'+Y','HorizontalAlignment','center','FontSize',9);
    text(1.05,0,'+X','HorizontalAlignment','left','FontSize',9);
    text(0,-1.055,'-Y','HorizontalAlignment','center','FontSize',9);
    text(-1.05,0,'-X','HorizontalAlignment','right','FontSize',9);
    legend([ha hb],{'alpha [0001] (measured)', ...
        'nearest beta <110> (reconstructed)'},'Location','southoutside');
    title(sprintf('Matched-pixel pole overlay | %d of %d pairs | median %.2f deg', ...
        np,n,qp(1)));
    saveFigure(folder,'E1_BOR_Alpha0001_Beta110_PairedOverlay.png');
catch MEPF
    warning('BOR overlay figure 18 failed: %s',MEPF.message);
end

% Fig. 19: paired normal-angle misfit vs full Burgers-variant misfit.
try
    figure('Name','E2 BOR plane vs full fit','Color','white');
    subplot(1,2,1);
    histogram(minPlaneDeg(validPlane),'BinWidth',0.5);
    grid on; xlabel('alpha c-axis vs closest beta {110} [deg]');
    ylabel('Matched EBSD pixels'); title('Necessary plane-normal match');
    subplot(1,2,2);
    histogram(fullBORfit_deg(validFull),'BinWidth',0.5);
    grid on; xlabel('Complete Burgers variant fit [deg]');
    ylabel('Matched EBSD pixels'); title('Full BOR (MTEX fit)');
    saveFigure(folder,'E2_BOR_PlaneNormal_vs_FullFit_Histograms.png');
catch MEhist
    warning('BOR misfit histograms failed: %s',MEhist.message);
end

% Fig. 20: a spatial map of basal normal -> nearest beta {110} deviation.
try
    figure('Name','E3 BOR plane misfit map','Color','white');
    plot(A,minPlaneDeg); set(gca,'Color',[0.12 0.12 0.12]);
    mtexColorbar('title','deg');
    title(sprintf('BOR basal-plane normal deviation | %.1f%% alpha coverage',coveragePct));
    saveFigure(folder,'E3_BOR_PlaneNormal_Misfit_Map.png');
catch MEmap
    warning('BOR plane-normal misfit map failed: %s',MEmap.message);
end

% Fig. 21: compare both diagnostics directly on identical matched pixels.
try
    figure('Name','E4 BOR metric comparison','Color','white');
    pp = find(pairedValid);
    if numel(pp)>25000
        pp = pp(round(linspace(1,numel(pp),25000)));
    end
    scatter(minPlaneDeg(pp),fullBORfit_deg(pp),7,'filled');
    grid on; xlabel('Basal vs beta-{110} normal deviation [deg]');
    ylabel('Full BOR variant misfit [deg]');
    title(sprintf('Paired fit comparison, %d pixels',numel(pp)));
    saveFigure(folder,'E4_BOR_Plane_vs_FullFit_Scatter.png');
catch MEscat
    warning('BOR misfit scatter failed: %s',MEscat.message);
end
end

function [xp,yp] = borLambertUpper(v)
% Unit-disk Lambert equal-area pole figure, upper hemisphere, specimen +X
% to the right and specimen +Y upwards. Poles and antipodes are the same.
x = v.x(:); y = v.y(:); z = v.z(:);
normXYZ = sqrt(x.^2+y.^2+z.^2);
x = x ./ normXYZ; y = y ./ normXYZ; z = z ./ normXYZ;
s = ones(numel(z),1); s(z<0) = -1;
x = x .* s; y = y .* s; z = z .* s;
den = sqrt(max(1+z,eps));
xp = x ./ den; yp = y ./ den;
end

%% ---------- SCHMID HELPERS -------------------------------------------------
function D = schmidAlphaAnalysis(alphaData,CSalpha,R,opt,folder)
% A force DIRECTION alone does not specify stress; SF is dimensionless.
% Data are computed for each indexed alpha pixel, not only BOR-matched pixels.
xyz = opt.loadAxisXYZ(:).';
if numel(xyz) ~= 3 || ~all(isfinite(xyz)) || norm(xyz) == 0
    error('SFopt.loadAxisXYZ must be a non-zero 3-element vector.');
end
xyz = xyz/norm(xyz);
loadDir = vector3d(xyz(1),xyz(2),xyz(3));

% MTEX documentation: sS.SchmidFactor(inv(orientation) * specimenDirection)
% makes the direction and slip systems share the CRYSTAL reference frame.
% antipodal symmetrization + |m| removes the arbitrary slip-sense sign.
ori = alphaData.orientations;
loadInCrystal = inv(ori) * loadDir;
sBas = slipSystem.basal(CSalpha).symmetrise('antipodal');
sPri = slipSystem.prismaticA(CSalpha).symmetrise('antipodal');
sPyr = slipSystem.pyramidalCA(CSalpha).symmetrise('antipodal');

% Reduce one family at a time to minimize memory usage for fine-step maps.
mBas = max(abs(sBas.SchmidFactor(loadInCrystal)),[],2);
mPri = max(abs(sPri.SchmidFactor(loadInCrystal)),[],2);
mPyr = max(abs(sPyr.SchmidFactor(loadInCrystal)),[],2);
mBas = mBas(:); mPri = mPri(:); mPyr = mPyr(:);

% |c-axis| vs load axis: 0 deg = c parallel to load, 90 deg = perpendicular.
cDirCrystal = Miller(0,0,0,1,CSalpha,'UVTW');
cAxisSpec = ori * cDirCrystal;
theta = angle(cAxisSpec,loadDir)./degree;
theta = theta(:);
theta = min(theta,180-theta);  % antipodal crystallographic axis
basalNormalRatio = cosd(theta).^2;  % sigma_n(basal)/sigma for tensile uniaxial load
mEasyA = max(mBas,mPri);            % geometry only, NOT an active slip prediction
matched = ismember(alphaData.id,R.rec.id);
matched = matched(:);

fprintf('\n--- ALPHA SCHMID FACTOR (-X external-load scenario) ---\n');
fprintf('Load axis in EBSD frame: [%g %g %g], direction towards LEFT if [-1 0 0]\n',xyz);
fprintf('Analyzed ALL indexed alpha: %d px | beta-BOR-matched alpha: %d (%.2f%%)\n', ...
    numel(mBas),nnz(matched),100*mean(matched));
fprintf('Uniaxial TENSION approximation. -X versus +X gives same absolute SF.\n');
fprintf('This computes potential slip geometry; no strain, force or local stress solved.\n');

% Quantify selection bias: beta-parent reconstruction generally covers a
% non-random subset of the alpha pixels, especially for coarse EBSD scans.
labels = {'ALL measured alpha','BOR-matched measured alpha','Alpha NOT BOR-matched'};
masks = {true(numel(mBas),1),matched,~matched};
% Optional user-marked candidate region: a geometric ROI is not itself proof
% that it meets crystallographic MTR segmentation criteria.
if ~isempty(opt.candidateMTRroiRect)
    r = opt.candidateMTRroiRect;
    if numel(r)~=4 || r(3)<=0 || r(4)<=0
        error('SFopt.candidateMTRroiRect needs [x0 y0 width height], with positive size.');
    end
    inROI = alphaData.x(:)>=r(1) & alphaData.x(:)<=r(1)+r(3) & ...
        alphaData.y(:)>=r(2) & alphaData.y(:)<=r(2)+r(4);
    labels = [labels, {'Candidate MTR ROI (user-defined)','Outside candidate MTR ROI'}];
    masks = [masks, {inROI,~inROI}];
    fprintf('Candidate ROI contains %d indexed alpha pixels. ROI is not automatic MTR detection.\n',nnz(inROI));
end
% MTR / macrozone pixels from RUN_MACROZONES (same .ctf), if available.
MZ = [];
try
    MZ = linkMacrozoneResults(alphaData,opt,folder);
catch MEmz
    warning('Macrozone link failed: %s',MEmz.message);
end
if ~isempty(MZ)
    labels = [labels, {'MTR pixels (macrozone tool)', ...
        sprintf('Macrozones >= %g um2',MZ.areaThreshold),'Alpha outside MTR'}];
    masks = [masks, {MZ.inMTR, MZ.id>0, ~MZ.inMTR}];
end
rows = cell(numel(labels),1);
for j=1:numel(labels)
    rows{j} = schmidSummaryRow(labels{j},masks{j},mBas,mPri,mPyr,mEasyA,theta);
end
summary = vertcat(rows{:});
writetable(summary,fullfile(folder,'F00_Schmid_Alpha_Summary.csv'));
disp(summary);

% Save numeric arrays in MATLAB format; use the CSV switch only if needed.
D = struct();
D.x_um = alphaData.x(:); D.y_um = alphaData.y(:);
D.id = alphaData.id(:); D.BORMatched = matched;
D.SF_Basal = mBas; D.SF_PrismaticA = mPri; D.SF_PyramidalCA = mPyr;
D.SF_MaxBasalOrPrismatic = mEasyA;
D.CAxisToLoad_deg = theta; D.BasalNormalStressOverSigma = basalNormalRatio;
D.loadAxisXYZ = xyz; D.summary = summary;

% Maps: whole measured alpha, not biased by whether parent beta was found.
plotAlphaSchmidMap(alphaData,mBas,[0 0.5],'m (absolute)', ...
    'F01 alpha: basal <a> Schmid factor','F01_SF_Basal_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mPri,[0 0.5],'m (absolute)', ...
    'F02 alpha: prismatic <a> Schmid factor','F02_SF_PrismaticA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mPyr,[0 0.5],'m (absolute)', ...
    'F03 alpha: pyramidal <c+a> Schmid factor','F03_SF_PyramidalCA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mEasyA,[0 0.5],'m (absolute)', ...
    'F04 alpha: max(basal, prismatic) | geometry only','F04_SF_MaxA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,theta,[0 90],'angle [deg]', ...
    'F05 alpha: acute angle between c-axis and load','F05_CAxis_Angle_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,basalNormalRatio,[0 1],'sigma_n / sigma', ...
    'F06 alpha: basal-plane normal stress / tensile sigma', ...
    'F06_Basal_NormalStressRatio_Alpha.png',xyz,folder);

% This is a plot only of alpha pixels that have reconstructed parent-beta
% orientations. It is a SUBSET, not the whole-sample alpha texture.
if any(matched)
    plotAlphaSchmidMap(alphaData(matched),mEasyA(matched),[0 0.5], ...
        'm (absolute)', ...
        sprintf('F07 BOR-matched alpha ONLY (%.1f%% of indexed alpha)',100*mean(matched)), ...
        'F07_SF_MaxA_BORmatched_Alpha_ONLY.png',xyz,folder);
end

% ---- Schmid factors PER MACROZONE (needs the macrozone link above) ----
if ~isempty(MZ)
    D.InMTR = MZ.inMTR; D.MacrozoneID = MZ.id; D.MacrozoneResultsFile = MZ.file;
    try
        Tmz = schmidPerMacrozone(MZ,mBas,mPri,mPyr,mEasyA,theta,basalNormalRatio,opt);
        writetable(Tmz,fullfile(folder,'F00b_Schmid_per_Macrozone.csv'));
        D.perMacrozone = Tmz;
        fprintf('\n--- SCHMID FACTOR PER MACROZONE (%d macrozones, file F00b) ---\n',height(Tmz));
        disp(Tmz(:,{'Macrozone','Area_um2','Basal_mean','Prismatic_mean', ...
            'MaxA_mean','CAxisToLoad_mean_deg','CAxisClass'}));
        cls = categorical(Tmz.CAxisClass);
        disp(table(categories(cls),countcats(cls),'VariableNames',{'CAxisClass','Macrozones'}));
    catch MEt
        warning('Per-macrozone Schmid table failed: %s',MEt.message);
    end
    inMZ = MZ.id > 0;
    if any(inMZ)
        plotAlphaSchmidMap(alphaData(inMZ),theta(inMZ),[0 90],'angle [deg]', ...
            sprintf('F14 macrozones >= %g um2: c-axis to load angle',MZ.areaThreshold), ...
            'F14_CAxis_Angle_Macrozones.png',xyz,folder);
        plotAlphaSchmidMap(alphaData(inMZ),mEasyA(inMZ),[0 0.5],'m (absolute)', ...
            sprintf('F15 macrozones >= %g um2: max(basal, prismatic) SF',MZ.areaThreshold), ...
            'F15_SF_MaxA_Macrozones.png',xyz,folder);
    end
    % F16: MTR in isolation vs the rest of the alpha phase
    try
        plotMTRvsNonMTR(MZ,mBas,mPri,theta,folder);
    catch MEh
        warning('MTR vs non-MTR plot failed: %s',MEh.message);
    end
    % F17: how each macrozone responds to the in-plane load direction
    if isfield(opt,'loadSweep') && opt.loadSweep && any(inMZ)
        try
            D.loadSweep = loadSweepMacrozones(ori,sBas,sPri,MZ,xyz,opt,folder,CSalpha);
        catch MEs
            warning('Load-direction sweep failed: %s',MEs.message);
        end
    end
end

% Comparable distribution curves for each slip family.
try
    e = linspace(0,0.5,51);
    center = 0.5*(e(1:end-1)+e(2:end));
    hb = histcounts(mBas,e,'Normalization','probability');
    hp = histcounts(mPri,e,'Normalization','probability');
    hc = histcounts(mPyr,e,'Normalization','probability');
    figure('Name','F08 Schmid distributions','Color','white');
    plot(center,hb,'LineWidth',1.7); hold on;
    plot(center,hp,'LineWidth',1.7);
    plot(center,hc,'LineWidth',1.7); hold off;
    xlabel('Maximum absolute Schmid factor per slip family');
    ylabel('Fraction of measured alpha pixels per bin');
    legend({'Basal <a>','Prismatic <a>','Pyramidal <c+a>'},'Location','best');
    xlim([0 0.5]); grid on;
    title('F08 All indexed alpha | geometric loading comparison');
    saveFigure(folder,'F08_SF_Histograms_Alpha.png');
catch ME
    warning('Schmid distribution plot failed: %s',ME.message);
end

% OPTIONAL: a tensile nominal stress lets us compute uniform-stress
% resolved shear tau = sigma*m. This is NOT a predicted local stress field.
if isfinite(opt.nominalStressMPa) && opt.nominalStressMPa > 0
    sigma = opt.nominalStressMPa;
    D.NominalStressMPa = sigma;
    D.TauBasal_MPa = sigma*mBas;
    D.TauPrismatic_MPa = sigma*mPri;
    D.TauPyramidalCA_MPa = sigma*mPyr;
    D.BasalNormalStress_MPa = sigma*basalNormalRatio;
    fprintf('Assumed UNIFORM tensile sigma = %.3f MPa (NOT derived from rotation speed).\n',sigma);
    plotAlphaSchmidMap(alphaData,D.TauBasal_MPa,[0 sigma/2],'|tau| [MPa]', ...
        'F09 basal |resolved shear| from assumed sigma','F09_Tau_Basal_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.TauPrismatic_MPa,[0 sigma/2],'|tau| [MPa]', ...
        'F10 prismatic |resolved shear| from assumed sigma','F10_Tau_Prismatic_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.TauPyramidalCA_MPa,[0 sigma/2],'|tau| [MPa]', ...
        'F11 pyramidal |resolved shear| from assumed sigma','F11_Tau_PyramidalCA_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.BasalNormalStress_MPa,[0 sigma], ...
        'sigma_n [MPa]', ...
        'F12 basal-plane normal stress from assumed sigma', ...
        'F12_Basal_NormalStress_Alpha.png',xyz,folder);
    if numel(opt.crssMPa) == 3 && all(isfinite(opt.crssMPa)) && all(opt.crssMPa>0)
        rb = D.TauBasal_MPa/opt.crssMPa(1);
        rp = D.TauPrismatic_MPa/opt.crssMPa(2);
        rc = D.TauPyramidalCA_MPa/opt.crssMPa(3);
        D.TauOverCRSS_Max = max([rb,rp,rc],[],2);
        D.CRSSFAMILY_MPa = opt.crssMPa;
        plotAlphaSchmidMap(alphaData,D.TauOverCRSS_Max, ...
            [0 max(1,max(D.TauOverCRSS_Max))],'max |tau| / CRSS', ...
            'F13 potential slip criterion (uniform stress + specified CRSS)', ...
            'F13_TauOverCRSS_Alpha.png',xyz,folder);
        fprintf('Fraction of alpha pixels with max(|tau|/CRSS) >= 1: %.2f%% (idealized)\n', ...
            100*mean(D.TauOverCRSS_Max>=1));
    else
        fprintf('CRSS unspecified: no slip-activation or plastic-strain prediction.\n');
    end
else
    fprintf('No tensile stress MPa given: stress/shear maps omitted; SF maps are valid.\n');
    if any(isfinite(opt.crssMPa))
        warning('CRSS provided but nominalStressMPa is NaN: cannot calculate tau/CRSS.');
    end
end

save(fullfile(folder,'F00_Schmid_Alpha_Results.mat'),'D','opt','-v7.3');
if opt.exportPixelCSV
    T = table(D.id,D.x_um,D.y_um,D.BORMatched,D.SF_Basal, ...
        D.SF_PrismaticA,D.SF_PyramidalCA,D.SF_MaxBasalOrPrismatic, ...
        D.CAxisToLoad_deg,D.BasalNormalStressOverSigma, ...
        'VariableNames',{'EBSD_ID','X_um','Y_um','BOR_matched', ...
        'SF_basal','SF_prismaticA','SF_pyramidalCA','SF_max_a', ...
        'CAxisToLoad_deg','BasalNormalStressOverSigma'});
    writetable(T,fullfile(folder,'F00_Schmid_Alpha_PerPixel.csv'));
end
fprintf('Schmid analysis saved in: %s\n',folder);
end


%% ---------- LINK TO MACROZONE (MTR) RESULTS -------------------------------
function MZ = linkMacrozoneResults(alphaData,opt,folder)
% Maps every alpha pixel of this MTEX analysis onto the grid of the
% macrozone analysis (RUN_MACROZONES / Read_CTF_Macrozones) of the SAME .ctf.
% Both grids come from the .ctf X,Y coordinates; the mapping is verified by
% checking that the alpha pixels land on Ti-Hex pixels of the macrozone grid.
MZ = [];
src = '';
if isfield(opt,'macrozoneResults'), src = opt.macrozoneResults; end
if isempty(src), return; end
[ctfDir,ctfName] = fileparts(opt.ctfFile);
if strcmpi(src,'run')
    % ---- built-in MTR detection on the same .ctf ----
    fprintf('\n--- MTR DETECTION (Fernandez Silva et al. 2023) on %s ---\n',ctfName);
    G = mtr_ReadCTF(opt.ctfFile,'auto');
    Pm = opt.MTR;
    if ischar(Pm.gridsize)
        g = round(Pm.windowUm/G.step); if mod(g,2)==0, g = g+1; end
        Pm.gridsize = max(3,g);
    end
    fprintf('  window %d px (%.1f um), %g deg, %g %%, area >= %g um2\n', ...
        Pm.gridsize,Pm.gridsize*G.step,Pm.critDisor,100*Pm.critFraction,Pm.areaThreshold);
    M = mtr_DetectMacrozones(G,Pm);
    G.G = []; G.Q = [];                          % keep the saved file small
    f = fullfile(folder,'MTR_Results.mat');
    ebsd = G; R = M; %#ok<NASGU>
    save(f,'ebsd','R','-v7.3');
    writetable(M.DataTable,fullfile(folder,'MTR_Macrozones_data.csv'));
    writetable(M.Statistics,fullfile(folder,'MTR_Macrozones_statistics.csv'),'WriteRowNames',true);
    try
        figure('Name','MTR macrozones','Color','white');
        imshow(M.MacroNumberColor_count,[]);
        title(sprintf('MTR macrozones >= %g um^2 (built-in detection): %d | MTR %.1f %%', ...
            Pm.areaThreshold,numel(M.kept),M.MTR),'Interpreter','none');
        saveFigure(folder,'F16a_MTR_Macrozones_map.png');
    catch
    end
elseif strcmpi(src,'ask')
    [a,b] = uigetfile(fullfile(ctfDir,'*.mat'), ...
        'Select Results.mat from RUN_MACROZONES (same .ctf)');
    if isequal(a,0), return; end
    f = fullfile(b,a);
elseif strcmpi(src,'auto')
    L = dir(fullfile(ctfDir,[ctfName '_Macrozones_*'],'Results.mat'));
    if isempty(L)
        fprintf(['Macrozone link: no %s_Macrozones_*\\Results.mat next to the .ctf.\n' ...
            '  Run RUN_MACROZONES on this .ctf first, or set SFopt.macrozoneResults.\n'],ctfName);
        return
    end
    [~,k] = max([L.datenum]);
    f = fullfile(L(k).folder,L(k).name);
    if numel(L) > 1
        fprintf('Macrozone link: %d result folders found, using the newest: %s\n',numel(L),L(k).folder);
    end
else
    f = src;
end
if ~strcmpi(src,'run')
    S = load(f,'ebsd','R');
    G = S.ebsd; M = S.R;
end
gName = char(G.name(:).');
if ~strcmp(gName,ctfName)
    warning('Macrozone results come from "%s", but this analysis uses "%s".',gName,ctfName);
end

% grid origin from the same .ctf (raw X,Y), as in Read_CTF_Macrozones
[minX,minY] = ctfMinXY(opt.ctfFile);
x = alphaData.x(:); y = alphaData.y(:);
% candidates: coordinate sign (MTEX may store -x/-y) x grid flips
cand = {}; names = {};
for sx = [1 -1]
    for sy = [1 -1]
        c0 = round((sx*x - minX)/G.step)  + 1;
        r0 = round((sy*y - minY)/G.stepY) + 1;
        base = sprintf('x%s y%s',char('+'*(sx>0)+'-'*(sx<0)),char('+'*(sy>0)+'-'*(sy<0)));
        cand(end+1,:) = {r0, c0};                    names{end+1} = base; %#ok<AGROW>
        cand(end+1,:) = {G.ny+1-r0, c0};             names{end+1} = [base ', flip rows']; %#ok<AGROW>
        cand(end+1,:) = {r0, G.nx+1-c0};             names{end+1} = [base ', flip cols']; %#ok<AGROW>
        cand(end+1,:) = {G.ny+1-r0, G.nx+1-c0};      names{end+1} = [base ', flip both']; %#ok<AGROW>
    end
end
score = zeros(size(cand,1),1);
for k = 1:size(cand,1)
    rr = cand{k,1}; cc = cand{k,2};
    ok = rr>=1 & rr<=G.ny & cc>=1 & cc<=G.nx;
    if ~any(ok), continue; end
    score(k) = mean(G.indexed(sub2ind([G.ny G.nx],rr(ok),cc(ok)))) * mean(ok);
end
[best,k] = max(score);
fprintf('Macrozone link: %s | grid match %.1f%% (mapping: %s)\n',f,100*best,names{k});
if best < 0.95
    warning(['Macrozone link: only %.1f%% of alpha pixels fall on Ti-Hex pixels of the ' ...
        'macrozone grid. Not the same map/crop? MTR rows skipped.'],100*best);
    MZ = []; return
end
rr = cand{k,1}; cc = cand{k,2};
ok = rr>=1 & rr<=G.ny & cc>=1 & cc<=G.nx;
lin = sub2ind([G.ny G.nx],rr(ok),cc(ok));
n = numel(x);
MZ.inMTR = false(n,1); MZ.inMTR(ok) = M.ismacro(lin);
MZ.id = zeros(n,1);    MZ.id(ok) = M.MacroNumber_count(lin);
MZ.dataTable = M.DataTable;
MZ.areaThreshold = M.P.areaThreshold;
MZ.file = f; MZ.gridMatch = best; MZ.flip = names{k};
fprintf('  alpha pixels in MTR: %d (%.1f%%) | in macrozones >= %g um2: %d (%.1f%%)\n', ...
    nnz(MZ.inMTR),100*mean(MZ.inMTR),MZ.areaThreshold,nnz(MZ.id>0),100*mean(MZ.id>0));
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
figure('Name','F16 MTR vs non-MTR','Color','white','Position',[80 80 1350 420]);
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
sgtitle('F16 Schmid factor of MTRs in isolation vs the rest of alpha (same load axis)');
saveFigure(folder,'F16_SF_MTR_vs_nonMTR.png');
end

function S = loadSweepMacrozones(ori,sBas,sPri,MZ,xyz,opt,folder,CSalpha)
% F17: mean basal / prismatic SF and c-axis angle of every macrozone as the
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
figure('Name','F17 load-direction sweep','Color','white','Position',[80 80 1350 420]);
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
sgtitle('F17 Sensitivity of each macrozone to the in-plane load direction (uniaxial, geometric)');
saveFigure(folder,'F17_LoadSweep_Macrozones.png');
% table: one row per angle, columns per macrozone
T = table(ang(:),'VariableNames',{'LoadAngle_deg'});
for z = 1:nZ
    T.(sprintf('MZ%d_Basal',ids(z))) = Bm(:,z);
    T.(sprintf('MZ%d_Prism',ids(z))) = Pm(:,z);
    T.(sprintf('MZ%d_MaxA',ids(z)))  = Am(:,z);
    T.(sprintf('MZ%d_CAxis_deg',ids(z))) = Th(:,z);
end
writetable(T,fullfile(folder,'F00c_LoadSweep_Macrozones.csv'));
% per macrozone: softest / hardest in-plane direction
[mx,ix] = max(Am,[],1); [mn,in] = min(Am,[],1);
[~,i0] = min(abs(ang - a0));
R = table(ids,Am(i0,:).',mx.',ang(ix).',mn.',ang(in).',max(Am,[],1).'-min(Am,[],1).', ...
    'VariableNames',{'Macrozone','MaxA_at_actual_load','MaxA_max','Angle_of_max_deg', ...
    'MaxA_min','Angle_of_min_deg','MaxA_range'});
writetable(R,fullfile(folder,'F00d_LoadSensitivity_Macrozones.csv'));
fprintf('\n--- LOAD-DIRECTION SENSITIVITY PER MACROZONE (F17, F00c/F00d) ---\n');
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


%% ========== SIX MATCHED ALPHA / RECONSTRUCTED-BETA POLE FIGURES =========
function figH = plotSixMatchedPFs(oriAlpha,oriBeta,csAlpha,csBeta,halfwidth,coveragePct,nBetaFragments,outputFolder)
% Explicit six-axis rendering for MATLAB + MTEX.
% A MATLAB subplot per PF is more reliable than a single newMtexFigure
% populated with alternating nextAxis commands.
% Row 1: measured alpha (0002), (10-10), (11-20).
% Row 2: reconstructed beta (110), (200), (111).
% IMPORTANT: oriAlpha and oriBeta must come from the SAME matched EBSD IDs.

    if isempty(oriAlpha) || isempty(oriBeta)
        error('Figure D4: empty matched alpha / reconstructed-beta populations.');
    end

    ha = [Miller(0,0,0,2,csAlpha), ...
          Miller(1,0,-1,0,csAlpha), ...
          Miller(1,1,-2,0,csAlpha)];
    hb = [Miller(1,1,0,csBeta), Miller(2,0,0,csBeta), ...
          Miller(1,1,1,csBeta)];
    la = {'{0002}','{10-10}','{11-20}'};
    lb = {'{110}','{200}','{111}'};

    odfA = calcDensity(oriAlpha,'halfwidth',halfwidth,'silent');
    odfB = calcDensity(oriBeta, 'halfwidth',halfwidth,'silent');

    % Maxima are computed from the SAME ODFs used to plot these 6 PFs.
    % MTEX: pdf = calcPDF(odf,h), then [maxMRD, specimenDirection] = max(pdf).
    % Do NOT substitute max(odf) (3D orientation density) for a PF peak.
    maxA = nan(1,3); maxB = nan(1,3);
    dirA = nan(3,3); dirB = nan(3,3);
    tiltA = nan(1,3); tiltB = nan(1,3);
    azA = nan(1,3); azB = nan(1,3);
    for k = 1:3
        [maxA(k),dirA(k,:),tiltA(k),azA(k)] = pfMaxWithDirection(odfA,ha(k));
        [maxB(k),dirB(k,:),tiltB(k),azB(k)] = pfMaxWithDirection(odfB,hb(k));
    end

    % Preserve clear labels for each figure: exact family and phase.
    familyNames = [la(:); lb(:)];
    phaseNames = [repmat({'Measured alpha'},3,1); repmat({'Reconstructed beta'},3,1)];
    maxValues = [maxA(:); maxB(:)];
    xyz = [dirA;dirB];
    tiltValues = [tiltA(:);tiltB(:)];
    azValues = [azA(:);azB(:)];
    nMatched = [repmat(numel(oriAlpha),3,1);repmat(numel(oriBeta),3,1)];
    halfwidthDeg = repmat(halfwidth/degree,6,1);
    pfStats = table(phaseNames,familyNames,maxValues,xyz(:,1),xyz(:,2),xyz(:,3), ...
        tiltValues,azValues,nMatched,halfwidthDeg, ...
        'VariableNames',{'Phase','PoleFamily','MaxMRD', ...
        'PeakDir_X','PeakDir_Y','PeakDir_Z','TiltFromZ_deg', ...
        'AzimuthFromX_deg','MatchedPixels','ODFHalfwidth_deg'});
    fprintf('\n===== FIGURE D4: SIX PF MAXIMA (MRD) =====\n');
    disp(pfStats);
    fprintf(['Directions refer to the EBSD specimen reference frame.\n' ...
             'Tilt: 0 deg = +Z; azimuth: 0 deg = +X, 90 deg = +Y.\n' ...
             'For antipodal PFs, peaks are reported in the upper hemisphere.\n']);
    try
        writetable(pfStats,fullfile(outputFolder, ...
            'D4_PF_MaxMRD_and_PeakDirections.csv'));
    catch MEcsv
        warning('Could not write PF maximum CSV: %s',MEcsv.message);
    end
    % Force undocked windows: MATLAB R2025a+ defaults figures to tabs,
    % which can make Coverage Figure 13 appear 'on top of' PF Figure D4.
    % An explicit figure handle keeps all six PF axes bound to THIS window.
    figH = figure('Name','D4 | SIX alpha-beta matched pole figures', ...
        'NumberTitle','off','Color','white','Visible','on', ...
        'WindowStyle','normal','Tag','Ti64_SixPF_Fig12', ...
        'Units','pixels','Position',[50 50 1600 850]);
    axs = gobjects(2,3);

    % Each PDF explicitly targets its OWN MATLAB axes.
    for k=1:3
        figure(figH);
        axs(1,k) = subplot(2,3,k,'Parent',figH);
        plotPDF(odfA,ha(k),'antipodal','contourf','parent',axs(1,k));
        title(axs(1,k),sprintf('Measured alpha %s | max %.2f MRD',la{k},maxA(k)), ...
            'FontWeight','bold','Interpreter','none','FontSize',11);

        figure(figH);
        axs(2,k) = subplot(2,3,3+k,'Parent',figH);
        plotPDF(odfB,hb(k),'antipodal','contourf','parent',axs(2,k));
        title(axs(2,k),sprintf('Reconstructed beta %s | max %.2f MRD',lb{k},maxB(k)), ...
            'FontWeight','bold','Interpreter','none','FontSize',11);
    end
    drawnow;
    if any(~isgraphics(axs(:),'axes'))
        error('Figure D4: MTEX removed or replaced one or more subplot axes.');
    end

    % Phase-specific MRD scales: same scale for the three alpha plots,
    % another same scale for the three beta plots. Avoid independent auto CLim.
    alphaRanges = zeros(1,3);
    betaRanges  = zeros(1,3);
    for k=1:3
        localA = get(axs(1,k),'CLim');
        localB = get(axs(2,k),'CLim');
        alphaRanges(k) = localA(2);
        betaRanges(k)  = localB(2);
    end
    aa=alphaRanges(isfinite(alphaRanges)); pa=maxA(isfinite(maxA));
    limAlpha = max([aa(:);pa(:)]);
    bb=betaRanges(isfinite(betaRanges)); pb=maxB(isfinite(maxB));
    limBeta  = max([bb(:);pb(:)]);
    if isempty(limAlpha) || limAlpha<=0, limAlpha=1; end
    if isempty(limBeta) || limBeta<=0, limBeta=1; end
    for k=1:3
        set(axs(1,k),'CLim',[0 limAlpha]);
        set(axs(2,k),'CLim',[0 limBeta]);
    end

    % One colorbar for each phase, on the right of its corresponding row.
    cbAlpha = colorbar(axs(1,3),'eastoutside');
    ylabel(cbAlpha,'MRD alpha');
    cbBeta = colorbar(axs(2,3),'eastoutside');
    ylabel(cbBeta,'MRD beta');
    try
        sgtitle(figH,sprintf(['alpha N = %d | beta N = %d | ', ...
            'alpha reconstruction coverage: %.1f%% | beta fragments: %d | ', ...
            'MTEX halfwidth %.1f deg (= Channel 5 %.0f deg)'], ...
            numel(oriAlpha),numel(oriBeta),coveragePct,nBetaFragments,halfwidth/degree,ch5EquivDeg(halfwidth)), ...
            'Interpreter','none','FontWeight','bold');
    catch
        % sgtitle may not be available in older MATLAB versions.
    end
    drawnow;

    outPng = fullfile(outputFolder,'D4_PF_Alpha_vs_Beta.png');
    outFig = fullfile(outputFolder,'D4_PF_Alpha_vs_Beta.fig');
    exportgraphics(figH,outPng,'Resolution',300);
    savefig(figH,outFig);
    fprintf('FIGURE D4 OK: SIX independent PFs displayed.\n');
    fprintf('Figure D4 PNG: %s\nFigure D4 FIG: %s\n',outPng,outFig);
    figure(figH); drawnow;
end

%% -------- Maximum of an MTEX pole density function ----------------------
function [peakMRD,peakXYZ,tiltDeg,azDeg] = pfMaxWithDirection(odf,h)
% Evaluates the maximum of the spherical POLE DENSITY FUNCTION, not the ODF.
% [value,direction] = max(calcPDF(odf,h)) is supported by MTEX S2Fun.max.
% The corresponding antipodal vector is flipped into the upper hemisphere,
% i.e., +Z-facing direction in the EBSD specimen reference frame.
    peakMRD = NaN;
    peakXYZ = [NaN NaN NaN];
    tiltDeg = NaN;
    azDeg = NaN;
    try
        pdfFun = calcPDF(odf,h);
        [peakMRD,peakDir] = max(pdfFun);
        peakMRD = double(peakMRD(1));
        if ~isfinite(peakMRD), error('Nonfinite maximum MRD.'); end
        xyz = [double(peakDir.x),double(peakDir.y),double(peakDir.z)];
        xyz = xyz(1,:);
        if xyz(3)<0, xyz = -xyz; end % antipodal equivalent
        xyz = xyz/norm(xyz);
        peakXYZ = xyz;
        tiltDeg = acos(max(-1,min(1,xyz(3))))/degree;
        if hypot(xyz(1),xyz(2)) > 1e-10
            azDeg = mod(atan2(xyz(2),xyz(1))/degree,360);
        end
    catch ME
        warning('Could not calculate a PF maximum/direction: %s',ME.message);
    end
end


%% ===== NINE POLE FIGURES: measured alpha, measured beta, reconstructed beta =====
function figH = plotNinePF_MeasuredAndReconstructedBeta( ...
    oriAlphaMatched,oriBetaMeasured,oriBetaReconstructed, ...
    csAlpha,csBeta,halfwidth,coveragePct,outFolder)
% Figure D5 is independent of Figure D4 and Figure 13 (own handle + tag).
% Row 1: measured alpha from the MATCHED reconstructed-alpha pixels ONLY.
% Row 2: ALL beta pixels MEASURED in the original CTF (NOT reconstructed).
% Row 3: reconstructed parent beta from the SAME matched alpha pixels as row 1.
% WARNING: row 2 is a different spatial population from rows 1 and 3.
% The measured-beta pixels are spatially correlated; counting each as an
% independent grain overstates effective texture sample size.

    if isempty(oriAlphaMatched) || isempty(oriBetaReconstructed)
        error('Figure D5 needs both matched alpha and reconstructed beta.');
    end

    ha = [Miller(0,0,0,2,csAlpha), ...
          Miller(1,0,-1,0,csAlpha), ...
          Miller(1,1,-2,0,csAlpha)];
    hb = [Miller(1,1,0,csBeta),Miller(2,0,0,csBeta), ...
          Miller(1,1,1,csBeta)];
    hByRow = {ha,hb,hb};
    labelsByRow = {{'{0002}','{10-10}','{11-20}'}, ...
                   {'{110}','{200}','{111}'}, ...
                   {'{110}','{200}','{111}'}};
    rowNames = {'Measured alpha', ...
                'Measured beta', ...
                'Reconstructed parent beta'};
    oriByRow = {oriAlphaMatched,oriBetaMeasured,oriBetaReconstructed};
    pixelCounts = [numel(oriAlphaMatched),numel(oriBetaMeasured), ...
                   numel(oriBetaReconstructed)];
    odfByRow = cell(1,3);
    available = false(1,3);
    for r=1:3
        if pixelCounts(r)==0
            warning('Figure D5: %s has ZERO pixels in the CTF.',rowNames{r});
            continue;
        end
        try
            odfByRow{r} = calcDensity(oriByRow{r}, ...
                'halfwidth',halfwidth,'silent');
            available(r) = true;
        catch MEodf
            warning('Figure D5: ODF failed for %s: %s', ...
                rowNames{r},MEodf.message);
        end
    end
    if ~available(1) || ~available(3)
        error('Figure D5: failed to compute alpha or reconstructed-beta ODF.');
    end
    if pixelCounts(2)>0 && pixelCounts(2)<100
        warning(['Only %d measured-beta EBSD pixels. Measured-beta PF ', ...
                 'is not a robust texture estimate.'],pixelCounts(2));
    end

    figH = figure('Name','D5 | NINE PF - alpha / measured beta / reconstructed beta', ...
        'NumberTitle','off','Color','white','Visible','on', ...
        'WindowStyle','normal','Tag','Ti64_NinePF_Fig25', ...
        'Units','pixels','Position',[120 90 1580 1000]);
    axs = gobjects(3,3);
    phaseNames = cell(9,1); familyNames = cell(9,1);
    maxValues = nan(9,1); xyz = nan(9,3);
    tiltVals = nan(9,1); azVals = nan(9,1);
    pixVals = zeros(9,1); halfVals = repmat(halfwidth/degree,9,1);

    for r=1:3
        for k=1:3
            q=(r-1)*3+k;
            phaseNames{q}=rowNames{r};
            familyNames{q}=labelsByRow{r}{k};
            pixVals(q)=pixelCounts(r);
            figure(figH);
            axs(r,k) = subplot(3,3,q,'Parent',figH);
            if available(r)
                [maxValues(q),xyz(q,:),tiltVals(q),azVals(q)] = ...
                    pfMaxWithDirection(odfByRow{r},hByRow{r}(k));
                plotPDF(odfByRow{r},hByRow{r}(k), ...
                    'antipodal','contourf','parent',axs(r,k));
                title(axs(r,k),sprintf('%s %s | max %.2f MRD', ...
                    rowNames{r},familyNames{q},maxValues(q)), ...
                    'FontSize',9,'FontWeight','bold','Interpreter','none');
            else
                axis(axs(r,k),'off');
                text(axs(r,k),0.5,0.5,'No indexed measured beta', ...
                    'Units','normalized','HorizontalAlignment','center', ...
                    'FontSize',10);
                title(axs(r,k),sprintf('%s %s | no data', ...
                    rowNames{r},familyNames{q}), ...
                    'Interpreter','none','FontSize',9);
            end
        end
    end
    drawnow;
    if any(~isgraphics(axs(:),'axes'))
        error('Figure D5: plotPDF replaced one or more of the nine axes.');
    end

    % Consistent MRD color limit within each phase/source ROW, NOT across
    % different populations: absolute MRD amplitudes remain in plot titles.
    for r=1:3
        if ~available(r), continue; end
        localUpper = zeros(1,3);
        for k=1:3
            lim = get(axs(r,k),'CLim');
            localUpper(k)=lim(2);
        end
        p=maxValues((r-1)*3+(1:3));
        % MATLAB can return one limit array as a column and another as
        % a row.  Concatenate COLUMN vectors to prevent horzcat errors.
        validLimits = localUpper(isfinite(localUpper));
        validPeaks = p(isfinite(p));
        allCandidates = [validLimits(:);validPeaks(:)];
        if isempty(allCandidates)
            upper = 1;
        else
            upper = max(allCandidates);
        end
        if isempty(upper) || ~isfinite(upper) || upper<=0, upper=1; end
        for k=1:3
            set(axs(r,k),'CLim',[0 upper]);
        end
        cb=colorbar(axs(r,3),'eastoutside');
        ylabel(cb,sprintf('MRD (row %d)',r));
    end

    try
        sgtitle(figH,sprintf(['Alpha: %d px | beta reconstruction coverage %.2f%% | ', ...
            'Beta measured: %d pixels | Beta reconstructed: %d pixels | ', ...
            'halfwidth %.1f deg (= Channel 5 %.0f deg)'], ...
            pixelCounts(1),coveragePct,pixelCounts(2),pixelCounts(3), ...
            halfwidth/degree,ch5EquivDeg(halfwidth)), ...
            'Interpreter','none','FontWeight','bold');
    catch
    end
    drawnow;

    T=table(phaseNames,familyNames,maxValues, ...
        xyz(:,1),xyz(:,2),xyz(:,3),tiltVals,azVals,pixVals,halfVals, ...
        'VariableNames',{'PhaseSource','PoleFamily','MaxMRD', ...
        'PeakDir_X','PeakDir_Y','PeakDir_Z', ...
        'TiltFromZ_deg','AzimuthFromX_deg','EBSD_Pixels','ODFHalfwidth_deg'});
    fprintf('\n===== FIGURE D5: NINE PF MAXIMA (MRD) =====\n');
    disp(T);
    fprintf(['Row 2 uses ALL originally MEASURED beta pixels; rows 1 and 3 ', ...
             'use only corresponding reconstructed alpha pixels.\n']);
    try
        writetable(T,fullfile(outFolder,'D5_PF_MaxMRD_and_PeakDirections.csv'));
    catch MEcsv
        warning('Figure D5 CSV export failed: %s',MEcsv.message);
    end
    pngFile=fullfile(outFolder,'D5_PF_Alpha_BetaMeasured_BetaReconstructed.png');
    figFile=fullfile(outFolder,'D5_PF_Alpha_BetaMeasured_BetaReconstructed.fig');
    % Save a handle explicitly: never use gcf (the coverage map may be active).
    try
        exportgraphics(figH,pngFile,'Resolution',300);
    catch MEpng
        warning('Figure D5 PNG export failed: %s',MEpng.message);
    end
    try
        savefig(figH,figFile);
    catch MEfig
        warning('Figure D5 FIG export failed: %s',MEfig.message);
    end
    fprintf('Figure D5: %s\n',pngFile);
    fprintf('Figure D5 editable FIG: %s\n',figFile);
end


%% ===== STANDALONE ORIGINAL MEASURED BETA: 3 PF + 3 IPF =====
function [figPF,figIPF] = plotMeasuredBetaPFandIPF(betaMeasured,halfwidth,showPF,outFolder)
% betaMeasured MUST be EBSD('Titanium cubic') from the original CTF.
% Figure D2: pole figures of {110}, {200}, {111} beta.
% Figure D3: inverse pole figures for the specimen directions X, Y, Z.
% These describe the original retained beta, NOT reconstructed prior beta.
% All three panels in each figure have their own explicit axes.
    figPF=[]; figIPF=[];
    n=length(betaMeasured);
    if n==0
        warning('No measured beta is indexed in the original CTF: figures 26/27 skipped.');
        return;
    end
    cs=betaMeasured.CS;
    odf=calcDensity(betaMeasured.orientations,'halfwidth',halfwidth,'silent');
    h=[Miller(1,1,0,cs),Miller(2,0,0,cs),Miller(1,1,1,cs)];
    labels={'{110}','{200}','{111}'};
    peakVals=nan(3,1); xyz=nan(3,3);
    tilt=nan(3,1); az=nan(3,1);

    %% Figure D2 - measured beta ONLY: three pole figures (= row 2 of fig 25).
    if showPF
    figPF=figure('Name','D2 | ORIGINAL MEASURED beta - 3 PFs', ...
        'NumberTitle','off','Color','white','Visible','on', ...
        'WindowStyle','normal','Tag','Ti64_MeasuredBeta_3PF_Fig26', ...
        'Units','pixels','Position',[180 160 1450 560]);
    axPF=gobjects(1,3);
    for k=1:3
        axPF(k)=subplot(1,3,k,'Parent',figPF);
        [peakVals(k),xyz(k,:),tilt(k),az(k)]=pfMaxWithDirection(odf,h(k));
        plotPDF(odf,h(k),'antipodal','contourf','parent',axPF(k));
        title(axPF(k),sprintf('Measured beta %s | max %.2f MRD', ...
            labels{k},peakVals(k)),'FontSize',11,'Interpreter','none');
    end
    drawnow;
    if any(~isgraphics(axPF,'axes'))
        warning('Figure D2: unexpected axes replacement by MTEX.');
    end
    tops=nan(3,1);
    for k=1:3
        lim=get(axPF(k),'CLim');
        tops(k)=lim(2);
    end
    allM=[tops(isfinite(tops));peakVals(isfinite(peakVals))];
    if isempty(allM), scale=1; else, scale=max(allM); end
    if ~isfinite(scale) || scale<=0, scale=1; end
    for k=1:3, set(axPF(k),'CLim',[0 scale]); end
    c=colorbar(axPF(3),'eastoutside'); ylabel(c,'MRD');
    try
        sgtitle(figPF,sprintf('ORIGINAL measured beta (CTF): %d indexed pixels | PF halfwidth %.1f deg', ...
            n,halfwidth/degree),'FontWeight','bold');
    catch
    end
    drawnow;
    T=table(labels(:),peakVals,xyz(:,1),xyz(:,2),xyz(:,3),tilt,az, ...
        repmat(n,3,1),repmat(halfwidth/degree,3,1), ...
        'VariableNames',{'PoleFamily','MaxMRD','PeakDir_X','PeakDir_Y', ...
        'PeakDir_Z','TiltFromZ_deg','AzimuthFromX_deg', ...
        'EBSD_Pixels','ODFHalfwidth_deg'});
    disp('===== ORIGINAL MEASURED BETA - THREE POLE FIGURE PEAKS =====');
    disp(T);
    writetable(T,fullfile(outFolder,'D2_MeasuredBeta_PF_MaxMRD.csv'));
    exportgraphics(figPF,fullfile(outFolder,'D2_MeasuredBeta_3PF.png'),'Resolution',300);
    savefig(figPF,fullfile(outFolder,'D2_MeasuredBeta_3PF.fig'));
    end

    %% Figure D3 - measured beta ONLY: three INVERSE pole figures.
    figIPF=figure('Name','D3 | ORIGINAL MEASURED beta - IPF-X/Y/Z', ...
        'NumberTitle','off','Color','white','Visible','on', ...
        'WindowStyle','normal','Tag','Ti64_MeasuredBeta_3IPF_Fig27', ...
        'Units','pixels','Position',[220 205 1450 560]);
    axesIPF=gobjects(1,3);
    dirs=[xvector,yvector,zvector];
    dlabels={'X','Y','Z'};
    for k=1:3
        axesIPF(k)=subplot(1,3,k,'Parent',figIPF);
        plotIPDF(odf,dirs(k),'contourf','parent',axesIPF(k));
        title(axesIPF(k),sprintf('Measured beta | IPF-%s',dlabels{k}), ...
            'FontSize',11,'Interpreter','none');
    end
    drawnow;
    limits=nan(3,1);
    for k=1:3
        lim=get(axesIPF(k),'CLim');
        limits(k)=lim(2);
    end
    upper=limits(isfinite(limits));
    if isempty(upper), upper=1; else, upper=max(upper); end
    if ~isfinite(upper) || upper<=0, upper=1; end
    for k=1:3, set(axesIPF(k),'CLim',[0 upper]); end
    c=colorbar(axesIPF(3),'eastoutside'); ylabel(c,'MRD');
    try
        sgtitle(figIPF,sprintf('Inverse pole density: MEASURED beta only | %d indexed pixels',n), ...
            'FontWeight','bold');
    catch
    end
    drawnow;
    exportgraphics(figIPF,fullfile(outFolder,'D3_MeasuredBeta_IPF_XYZ.png'),'Resolution',300);
    savefig(figIPF,fullfile(outFolder,'D3_MeasuredBeta_IPF_XYZ.fig'));
    fprintf('Figure D3: 3 measured beta IPF saved.\n');
end


%% ===== v4.9 THREE INDEPENDENT PF POPULATIONS ===============================
function [figAll,figMatched,figRec] = plotThreeAlphaBetaPopulations( ...
    alphaAll,alphaMatched,betaRec,csAlpha,csBeta,halfwidth,coverage,showRedundant,outFolder)
% 28: alpha ALL indexed pixels (CHANNEL 5 whole-phase comparison)
% 29: alpha MATCHED to reconstructed beta only
% 30: beta RECONSTRUCTED at these same measured alpha positions
% Each figure uses separate axes, same within-population color scale.
% No original measured-beta orientations are mixed into fig30.

    figAll = []; figMatched = []; figRec = [];
    if isempty(alphaAll)
        error('No indexed alpha data are available.');
    end
    nAll = length(alphaAll);
    nMatch = length(alphaMatched);
    nRec = length(betaRec);
    if nMatch ~= nRec
        error('Matched alpha and reconstructed beta pixel counts differ.');
    end

    [figAll,Tall] = threePFOnePopulation(alphaAll.orientations,csAlpha, ...
        halfwidth,'Whole-map measured alpha',1,'D1_AlphaWholeMap_3PF',outFolder);
    writetable(Tall,fullfile(outFolder,'D1_AlphaWholeMap_PF_MaxMRD.csv'));
    if ~showRedundant || nMatch == 0
        return;   % 29/30 repeat fig 12 (same pixels, same kernel)
    end
    [figMatched,Tmatched] = threePFOnePopulation(alphaMatched.orientations,csAlpha, ...
        halfwidth,'Matched measured alpha',9,'D9_AlphaMatched_3PF',outFolder);
    [figRec,Trec] = threePFOnePopulation(betaRec.orientations,csBeta, ...
        halfwidth,'Reconstructed parent beta',10,'D10_BetaReconstructed_3PF',outFolder);

    % Combined table, suitable for Excel comparison and research reports.
    T = [Tall; Tmatched; Trec];
    writetable(T,fullfile(outFolder,'D1_D9_D10_ThreePopulations_PF_MaxMRD.csv'));
    fprintf('\n===== THREE PF POPULATIONS SUMMARY =====\n');
    fprintf('All measured alpha pixels:        %d\n',nAll);
    fprintf('Matched measured alpha pixels:    %d\n',nMatch);
    fprintf('Reconstructed beta orientations:  %d\n',nRec);
    fprintf('Alpha reconstruction coverage: %.2f %%\n',coverage);
    disp(T(:,{'Population','PoleFamily','MaxMRD','EBSD_Pixels'}));
    fprintf(['COMPARE FIGURE 28 with CHANNEL 5 alpha phase only, same ROI, ' ...
             'coordinate frame, kernel and weighting.\n']);
    fprintf('Figures 29 and 30 are paired subsets, NOT whole-sample texture.\n');
    % Do not force figure 28 to front: existing fig 12/26/27 windows remain intact.
end

function [figH,T] = threePFOnePopulation(ori,cs,halfwidth,label,figNumber,stem,outFolder)
    if isempty(ori), error('Empty orientation population: %s',label); end
    if figNumber == 10
        hh = [Miller(1,1,0,cs),Miller(2,0,0,cs),Miller(1,1,1,cs)];
        families = {'{110}';'{200}';'{111}'};
    else
        hh = [Miller(0,0,0,2,cs),Miller(1,0,-1,0,cs), ...
              Miller(1,1,-2,0,cs)];
        families = {'{0002}';'{10-10}';'{11-20}'};
    end
    odf = calcDensity(ori,'halfwidth',halfwidth,'silent');
    peak = nan(3,1); xyz = nan(3,3); tilt=nan(3,1); az=nan(3,1);
    for k=1:3
        [peak(k),xyz(k,:),tilt(k),az(k)] = pfMaxWithDirection(odf,hh(k));
    end
    figH = figure('Name',sprintf('D%d | %s - three PFs',figNumber,label), ...
        'NumberTitle','off','Color','white','WindowStyle','normal', ...
        'Visible','on','Tag',sprintf('Ti64_PF_Population_%02d',figNumber), ...
        'Units','pixels','Position',[90 90 1480 570]);
    ax=gobjects(1,3);
    for k=1:3
        ax(k)=subplot(1,3,k,'Parent',figH);
        plotPDF(odf,hh(k),'antipodal','contourf','parent',ax(k));
        title(ax(k),sprintf('%s | max %.2f MRD',families{k},peak(k)), ...
            'FontWeight','bold','Interpreter','none','FontSize',11);
    end
    drawnow;
    if any(~isgraphics(ax,'axes'))
        warning('Figure %d: some PF axes are not available.',figNumber);
    end
    climTop=nan(3,1);
    for k=1:3
        lim=get(ax(k),'CLim'); climTop(k)=lim(2);
    end
    cand=[climTop(isfinite(climTop));peak(isfinite(peak))];
    if isempty(cand), upper=1; else, upper=max(cand); end
    if ~isfinite(upper) || upper<=0, upper=1; end
    for k=1:3, set(ax(k),'CLim',[0 upper]); end
    cb=colorbar(ax(3),'eastoutside'); ylabel(cb,'MRD');
    try
        sgtitle(figH,sprintf('%s | N = %d indexed orientations | MTEX halfwidth %.1f deg (= Channel 5 %.0f deg)', ...
            label,numel(ori),halfwidth/degree,ch5EquivDeg(halfwidth)),'FontWeight','bold');
    catch
    end
    drawnow;
    population=repmat({label},3,1);
    T=table(population,families,peak,xyz(:,1),xyz(:,2),xyz(:,3), ...
        tilt,az,repmat(numel(ori),3,1),repmat(halfwidth/degree,3,1), ...
        repmat(ch5EquivDeg(halfwidth),3,1), ...
        'VariableNames',{'Population','PoleFamily','MaxMRD', ...
        'PeakDir_X','PeakDir_Y','PeakDir_Z','TiltFromZ_deg', ...
        'AzimuthFromX_deg','EBSD_Pixels','ODFHalfwidth_deg','Channel5_equiv_deg'});
    writetable(T,fullfile(outFolder,[stem '_MaxMRD.csv']));
    exportgraphics(figH,fullfile(outFolder,[stem '.png']),'Resolution',300);
    savefig(figH,fullfile(outFolder,[stem '.fig']));
    fprintf('Figure %d generated: %s\n',figNumber,stem);
end
