%% Ti-6Al-4V: prior-beta reconstruction, MTEX 6.1.1  (v4, BOR overlay)
% v4 additions:
%   - Figure 18: overlaid α [0001] and NEAREST reconstructed β <110> poles,
%     evaluated ONLY for spatially corresponding EBSD pixels (same unique ID).
%   - Figures 19/20/21: distribution, map and scatter for BOR angular misfit.
%   - Pixel-by-pixel basal normal / {110}β alignment is a NECESSARY, NOT
%     SUFFICIENT, criterion for the full Burgers orientation relationship.
%   - Full BOR validation uses parentEBSD.fit from MTEX as an additional
%     variant-resolved misfit (not inferred from pole positions alone).
%   - Fix an undefined SFopt.crssMPa that could terminate Schmid analysis.
% Original v3 changes vs v2 (driven by the results of the 0.15 um map):
%   - Fragment pixel counts now use measured pixels (grainSize), not area/step^2
%     (area includes interior holes: 40076 "px" by area vs 33595 reconstructed px).
%   - Field-of-view criterion now uses the 90th percentile of ALL beta ECDs,
%     because interior-only statistics are biased low (big grains touch the edge).
%   - New reliability criteria: fraction of grains touching edge/holes and
%     fraction of physically implausible (tiny) beta grains.
%   - Pole figures: restored the pixel-weighted and grain-area-weighted beta PFs
%     of the original script (figures 10 and 11), shortened panel titles
%     (they overlapped), and warn based on the number of independent grains.
%   - Coverage summary pooled by alpha grain size (<3, 3-8, >=8 um).
%   - Empty-histogram handling; figure 03 only when the BOR is refined.
%
% NOTE: edited without being able to run MATLAB/MTEX. The main path of v2 ran
% on your data; the sensitivity and convergence modules are still untested.
%
% A reconstructed beta orientation is an INFERENCE from the BOR; it does not
% prove that each alpha grain came from a single recoverable beta parent.

clear; close all; clc;

%% 0. Start MTEX
if exist('EBSD','class') ~= 8
    startupCandidates = {'C:\mtex-6.1.1\startup_mtex.m','C:\mtex-6.1.1\startup.m'};
    started = false;
    for k = 1:numel(startupCandidates)
        if exist(startupCandidates{k},'file') == 2
            run(startupCandidates{k}); started = true; break
        end
    end
    if ~started, error('MTEX not found. Edit startupCandidates.'); end
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

pfHalfwidth = 10*degree;

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

try
    figure('Name','01 Alpha IPF-Z','Color','white');
    plot(alphaData,alphaKey.orientation2color(alphaData.orientations));
    set(gca,'Color',[0.12 0.12 0.12]);
    title([sampleName ' | measured alpha IPF-Z'],'Interpreter','none');
    saveFigure(outputFolder,'01_Alpha_IPFz.png');

    if ~isempty(R.gbFit0)
        figure('Name','02 Initial BOR fit','Color','white'); histogram(R.gbFit0,50);
        xlabel('BOR misfit per boundary [deg]'); ylabel('Neighbour pairs');
        title('Ideal BOR: fit on alpha/alpha boundaries');
        saveFigure(outputFolder,'02_BORfit_Initial.png');
    end
    if P.refineBOR && ~isempty(R.gbFit)
        figure('Name','03 BOR fit used','Color','white'); histogram(R.gbFit,50);
        xlabel('BOR misfit per boundary [deg]'); ylabel('Neighbour pairs');
        title('BOR used in the reconstruction');
        saveFigure(outputFolder,'03_BORfit_Used.png');
    end

    figure('Name','04 Beta before merging','Color','white');
    plot(R.parentBefore,betaKey.orientation2color(R.parentBefore.meanOrientation));
    hold on; plot(R.parentBefore.boundary,'lineWidth',0.25,'lineColor','white'); hold off;
    title('Beta fragments | before mergeSimilar/mergeInclusions');
    saveFigure(outputFolder,'04_Beta_BEFORE_merge.png');

    figure('Name','05 Final beta IPF-Z','Color','white');
    plot(R.parentGrains,betaKey.orientation2color(R.parentGrains.meanOrientation));
    hold on; plot(R.parentGrains.boundary,'lineWidth',0.2,'lineColor','white'); hold off;
    title('Reconstructed beta (job.parentGrains) | IPF-Z');
    saveFigure(outputFolder,'05_ParentBeta_IPFz.png');

    figure('Name','06 Beta boundaries on alpha','Color','white');
    plot(alphaData,alphaKey.orientation2color(alphaData.orientations));
    hold on; plot(R.parentGrains.boundary,'lineWidth',0.1,'lineColor','black'); hold off;
    title('Measured alpha | inferred beta boundaries');
    saveFigure(outputFolder,'06_BetaBoundaries_on_Alpha.png');

    if R.nRec > 0
        figure('Name','07 Beta per pixel IPF-Z','Color','white');
        plot(R.rec,betaKey.orientation2color(R.rec.orientations));
        title('Reconstructed beta per pixel (converted alpha pixels only)');
        saveFigure(outputFolder,'07_Beta_IPFz_PIXELS.png');

        figure('Name','08 BOR misfit per pixel','Color','white');
        plot(R.rec,R.fit_deg); mtexColorbar; setColorRange([0 5]);
        hold on; plot(R.parentGrains.boundary,'lineWidth',0.5,'lineColor','white'); hold off;
        title('Measured alpha vs reconstructed beta misfit [deg]');
        saveFigure(outputFolder,'08_Beta_PixelFit.png');
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
    figure('Name','13 Coverage map','Color','white');
    plot(ebsdRaw,cat); colormap(gca,cmap); setColorRange([0.5 4.5]);
    hold on;
    hl = gobjects(4,1);
    for k = 1:4, hl(k) = patch(NaN,NaN,cmap(k,:)); end
    legend(hl,{'Not indexed','Alpha NOT reconstructed','Alpha reconstructed','Measured beta'}, ...
        'Location','southoutside','Orientation','horizontal');
    hold off;
    title(sprintf('Coverage: %.1f%% of indexed alpha',R.coverage_pct));
    saveFigure(outputFolder,'13_Coverage_Map.png');
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
        figure('Name','14 Fragmentation','Color','white');
        edgesPx = [0.5 1.5 2.5 3.5 5.5 10.5 20.5 50.5 100.5 500.5 Inf];
        cnt = histcounts(S.nPx,edgesPx);
        bar(cnt); set(gca,'XTickLabel',{'1','2','3','4-5','6-10','11-20','21-50','51-100','101-500','>500'});
        xlabel('Pixels per beta grain'); ylabel('Number of grains');
        title('Pixels per reconstructed beta grain');
        saveFigure(outputFolder,'14_Fragmentation.png');

        figure('Name','09 ECD histograms','Color','white');
        subplot(2,2,1); plotHist(S.ecd,'All');
        subplot(2,2,2); plotHist(S.ecd(~S.edge),'Interior');
        subplot(2,2,3); plotHist(S.ecd(S.edge),'Touching edge/holes');
        subplot(2,2,4); plotHist(S.ecd(S.nPx<=3),'Fragments <= 3 px');
        saveFigure(outputFolder,'09_Beta_ECD_Histograms.png');
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

%% 12. Pole figures
if R.nRec > 0
    try
        hA = [Miller(0,0,0,2,CSalpha),Miller(1,0,-1,0,CSalpha),Miller(1,1,-2,0,CSalpha)];
        hB = [Miller(1,1,0,CSbeta),Miller(2,0,0,CSbeta),Miller(2,1,1,CSbeta)];
        tA = {'{0002}','{10-10}','{11-20}'};  tB = {'{110}','{200}','{211}'};

        % 12.1 beta PF, pixel-weighted (original figure 10)
        odfBetaPixel = calcDensity(R.rec.orientations,'halfwidth',pfHalfwidth,'silent');
        figure('Name','10 Beta PF pixel-weighted','Color','white');
        plotPDF(odfBetaPixel,hB,'antipodal','contourf'); mtexColorbar('title','mrd');
        saveFigure(outputFolder,'10_ParentBeta_PF_PixelWeighted.png');

        % 12.2 beta PF, one orientation per grain weighted by area (original figure 11)
        if S.n > 0
            odfBetaArea = calcDensity(S.ori,'weights',S.area,'halfwidth',pfHalfwidth,'silent');
            figure('Name','11 Beta PF grain-area-weighted','Color','white');
            plotPDF(odfBetaArea,hB,'antipodal','contourf'); mtexColorbar('title','mrd');
            saveFigure(outputFolder,'11_ParentBeta_PF_GrainAreaWeighted.png');
        end

        % 12.3 matched alpha vs beta (same pixels)
        odfAlphaMatched = calcDensity(alphaMatched.orientations,'halfwidth',pfHalfwidth,'silent');
        odfBetaMatched  = odfBetaPixel;
        figComp = newMtexFigure('layout',[3,2],'figSize','large');
        for iPF = 1:3
            plotPDF(odfAlphaMatched,hA(iPF),'antipodal','contourf');
            mtexTitle(sprintf('alpha %s | %d px',tA{iPF},R.nRec));
            nextAxis;
            plotPDF(odfBetaMatched,hB(iPF),'antipodal','contourf');
            mtexTitle(sprintf('beta %s | %d grains',tB{iPF},S.n));
            if iPF < 3, nextAxis; end
        end
        setColorRange('equal'); mtexColorbar('title','mrd'); drawNow(figComp);
        try, sgtitle(sprintf('Matched pixels: %.1f%% of the indexed alpha', ...
                R.coverage_pct)); catch, end
        saveFigure(outputFolder,'12_PF_Alpha_vs_Beta_MatchedPixels.png');

        % 12.4 reference: alpha of the WHOLE map (selection bias)
        odfAlphaAll = calcDensity(alphaData.orientations,'halfwidth',pfHalfwidth,'silent');
        figAll = newMtexFigure('layout',[1,3],'figSize','large');
        for iPF = 1:3
            plotPDF(odfAlphaAll,hA(iPF),'antipodal','contourf');
            mtexTitle(sprintf('alpha %s | all px',tA{iPF}));
            if iPF < 3, nextAxis; end
        end
        setColorRange('equal'); mtexColorbar('title','mrd'); drawNow(figAll);
        try, sgtitle(sprintf('Whole map: N = %d alpha pixels',length(alphaData))); catch, end
        saveFigure(outputFolder,'12b_PF_Alpha_AllPixels_Reference.png');
    catch ME
        warning('Pole figures not generated: %s',ME.message);
    end
    if R.coverage_pct < 50
        warning(['The beta PF describes only %.1f%% of the indexed alpha: it is NOT ', ...
            'the beta texture of the whole sample.'],R.coverage_pct);
    end
    if S.n < 100
        warning(['The beta PF is built from only %d independent orientations (grains). ', ...
            'Pixel counts do not add statistical weight: treat it as a list of ', ...
            'orientations, not as a texture estimate.'],S.n);
    end
end


%% 12b. PAIRED BOR VALIDATION: ALPHA (0001) vs NEAREST BETA {110}
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

%% 12c. SCHMID FACTOR: external load pointing LEFT (-X on EBSD map)
% Calculates SF for ALL measured alpha pixels, independently of BOR coverage.
% Reconstructed beta is a HISTORICAL parent, not today's alpha slip carrier.
% The Schmid factor predicts geometric resolved shear, NOT observed slip or
% a mechanical response. For genuine centrifugal loading, the full local
% stress tensor depends on specimen geometry, rotation and boundary conditions.
if runSchmidAnalysis
    try
        SFresults = schmidAlphaAnalysis(alphaData,CSalpha,R,SFopt,outputFolder);
    catch ME
        warning('Schmid-factor analysis failed: %s',ME.message);
        fprintf('Schmid detail: %s\n',getReport(ME,'extended','hyperlinks','off'));
    end
end

%% 13. Comparison with your original parameters
cmpRows = metricsRow(R,P,stepUm,'main');
if runBaselineComparison
    Rb = runReconstruction(ebsdRaw,alphaPhase,betaPhase,Pbase,stepUm);
    cmpRows = [cmpRows; metricsRow(Rb,Pbase,stepUm,'baseline_original')];
    fprintf('\n--- COMPARISON: new vs original parameters ---\n');
    disp(cmpRows);
    writetable(cmpRows,fullfile(outputFolder,'Comparison_main_vs_baseline.csv'));
end

%% 14. Parameter sensitivity (optional)
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

%% 15. Convergence versus step size (optional)
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
        figure('Name','16 Convergence','Color','white');
        loglog(Tconv.Step_um,Tconv.MedianECD_um,'o-',Tconv.Step_um,Tconv.MeanECD_um,'s-'); hold on;
        loglog(Tconv.Step_um,Tconv.Step_um*Tconv.MedianECD_um(1)/Tconv.Step_um(1),'k--');
        hold off; grid on;
        xlabel('Step size [um]'); ylabel('Beta ECD [um]');
        legend('Median','Mean','Slope 1 (size = k x pixel)','Location','northwest');
        title('If it follows the dashed line, the size is a pixel artefact');
        saveFigure(outputFolder,'16_Convergence.png');
    catch ME
        warning('Convergence figure not generated: %s',ME.message);
    end
end

%% 16. Save
settings = struct('P',P,'Pbase',Pbase,'Q',Q,'inputFile',fname,'stepUm',stepUm, ...
    'roiRect',roiRect,'downsampleFactor',downsampleFactor,'reliability',V, ...
    'runSchmidAnalysis',runSchmidAnalysis,'SFopt',SFopt, ...
    'runBORPairwiseCheck',runBORPairwiseCheck,'BORopt',BORopt);
Rsave = rmfield(R,intersect(fieldnames(R),{'job'}));
save(fullfile(outputFolder,'ParentBeta_v4_BOR_Schmid_Workspace.mat'),'Rsave','settings','-v7.3');
fprintf('\nDONE. Results in:\n%s\n',outputFolder);

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
writetable(Tsum,fullfile(folder,'18_BOR_Pairwise_Summary.csv'));
save(fullfile(folder,'18_BOR_Pairwise_AllPixels.mat'),'B','-v7.3');
if opt.exportPixelCSV
    Tout = table(idA,A.x(:),A.y(:),minPlaneDeg,fullBORfit_deg,double(nearestPole), ...
        'VariableNames',{'EBSD_ID','X_um','Y_um','PlaneNormalDeviation_deg', ...
        'FullBORVariantFit_deg','SelectedBeta110FamilyIndex'});
    writetable(Tout,fullfile(folder,'18_BOR_Pairwise_AllPixels.csv'));
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
    figure('Name','18 Alpha basal / nearest beta110 PF overlay','Color','white');
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
    saveFigure(folder,'18_BOR_Alpha0001_Beta110_PairedOverlay.png');
catch MEPF
    warning('BOR overlay figure 18 failed: %s',MEPF.message);
end

% Fig. 19: paired normal-angle misfit vs full Burgers-variant misfit.
try
    figure('Name','19 BOR plane vs full fit','Color','white');
    subplot(1,2,1);
    histogram(minPlaneDeg(validPlane),'BinWidth',0.5);
    grid on; xlabel('alpha c-axis vs closest beta {110} [deg]');
    ylabel('Matched EBSD pixels'); title('Necessary plane-normal match');
    subplot(1,2,2);
    histogram(fullBORfit_deg(validFull),'BinWidth',0.5);
    grid on; xlabel('Complete Burgers variant fit [deg]');
    ylabel('Matched EBSD pixels'); title('Full BOR (MTEX fit)');
    saveFigure(folder,'19_BOR_PlaneNormal_vs_FullFit_Histograms.png');
catch MEhist
    warning('BOR misfit histograms failed: %s',MEhist.message);
end

% Fig. 20: a spatial map of basal normal -> nearest beta {110} deviation.
try
    figure('Name','20 BOR plane misfit map','Color','white');
    plot(A,minPlaneDeg); set(gca,'Color',[0.12 0.12 0.12]);
    mtexColorbar('title','deg');
    title(sprintf('BOR basal-plane normal deviation | %.1f%% alpha coverage',coveragePct));
    saveFigure(folder,'20_BOR_PlaneNormal_Misfit_Map.png');
catch MEmap
    warning('BOR plane-normal misfit map failed: %s',MEmap.message);
end

% Fig. 21: compare both diagnostics directly on identical matched pixels.
try
    figure('Name','21 BOR metric comparison','Color','white');
    pp = find(pairedValid);
    if numel(pp)>25000
        pp = pp(round(linspace(1,numel(pp),25000)));
    end
    scatter(minPlaneDeg(pp),fullBORfit_deg(pp),7,'filled');
    grid on; xlabel('Basal vs beta-{110} normal deviation [deg]');
    ylabel('Full BOR variant misfit [deg]');
    title(sprintf('Paired fit comparison, %d pixels',numel(pp)));
    saveFigure(folder,'21_BOR_Plane_vs_FullFit_Scatter.png');
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
rows = cell(numel(labels),1);
for j=1:numel(labels)
    rows{j} = schmidSummaryRow(labels{j},masks{j},mBas,mPri,mPyr,mEasyA,theta);
end
summary = vertcat(rows{:});
writetable(summary,fullfile(folder,'17_Schmid_Alpha_Summary.csv'));
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
    '17a alpha: basal <a> Schmid factor','17a_SF_Basal_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mPri,[0 0.5],'m (absolute)', ...
    '17b alpha: prismatic <a> Schmid factor','17b_SF_PrismaticA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mPyr,[0 0.5],'m (absolute)', ...
    '17c alpha: pyramidal <c+a> Schmid factor','17c_SF_PyramidalCA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,mEasyA,[0 0.5],'m (absolute)', ...
    '17d alpha: max(basal, prismatic) | geometry only','17d_SF_MaxA_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,theta,[0 90],'angle [deg]', ...
    '17e alpha: acute angle between c-axis and load','17e_CAxis_Angle_Alpha.png',xyz,folder);
plotAlphaSchmidMap(alphaData,basalNormalRatio,[0 1],'sigma_n / sigma', ...
    '17f alpha: basal-plane normal stress / tensile sigma', ...
    '17f_Basal_NormalStressRatio_Alpha.png',xyz,folder);

% This is a plot only of alpha pixels that have reconstructed parent-beta
% orientations. It is a SUBSET, not the whole-sample alpha texture.
if any(matched)
    plotAlphaSchmidMap(alphaData(matched),mEasyA(matched),[0 0.5], ...
        'm (absolute)', ...
        sprintf('17g BOR-matched alpha ONLY (%.1f%% of indexed alpha)',100*mean(matched)), ...
        '17g_SF_MaxA_BORmatched_Alpha_ONLY.png',xyz,folder);
end

% Comparable distribution curves for each slip family.
try
    e = linspace(0,0.5,51);
    center = 0.5*(e(1:end-1)+e(2:end));
    hb = histcounts(mBas,e,'Normalization','probability');
    hp = histcounts(mPri,e,'Normalization','probability');
    hc = histcounts(mPyr,e,'Normalization','probability');
    figure('Name','17h Schmid distributions','Color','white');
    plot(center,hb,'LineWidth',1.7); hold on;
    plot(center,hp,'LineWidth',1.7);
    plot(center,hc,'LineWidth',1.7); hold off;
    xlabel('Maximum absolute Schmid factor per slip family');
    ylabel('Fraction of measured alpha pixels per bin');
    legend({'Basal <a>','Prismatic <a>','Pyramidal <c+a>'},'Location','best');
    xlim([0 0.5]); grid on;
    title('17h All indexed alpha | geometric loading comparison');
    saveFigure(folder,'17h_SF_Histograms_Alpha.png');
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
        '17i basal |resolved shear| from assumed sigma','17i_Tau_Basal_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.TauPrismatic_MPa,[0 sigma/2],'|tau| [MPa]', ...
        '17j prismatic |resolved shear| from assumed sigma','17j_Tau_Prismatic_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.TauPyramidalCA_MPa,[0 sigma/2],'|tau| [MPa]', ...
        '17k pyramidal |resolved shear| from assumed sigma','17k_Tau_PyramidalCA_Alpha.png',xyz,folder);
    plotAlphaSchmidMap(alphaData,D.BasalNormalStress_MPa,[0 sigma], ...
        'sigma_n [MPa]', ...
        '17l basal-plane normal stress from assumed sigma', ...
        '17l_Basal_NormalStress_Alpha.png',xyz,folder);
    if numel(opt.crssMPa) == 3 && all(isfinite(opt.crssMPa)) && all(opt.crssMPa>0)
        rb = D.TauBasal_MPa/opt.crssMPa(1);
        rp = D.TauPrismatic_MPa/opt.crssMPa(2);
        rc = D.TauPyramidalCA_MPa/opt.crssMPa(3);
        D.TauOverCRSS_Max = max([rb,rp,rc],[],2);
        D.CRSSFAMILY_MPa = opt.crssMPa;
        plotAlphaSchmidMap(alphaData,D.TauOverCRSS_Max, ...
            [0 max(1,max(D.TauOverCRSS_Max))],'max |tau| / CRSS', ...
            '17m potential slip criterion (uniform stress + specified CRSS)', ...
            '17m_TauOverCRSS_Alpha.png',xyz,folder);
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

save(fullfile(folder,'17_Schmid_Alpha_Results.mat'),'D','opt','-v7.3');
if opt.exportPixelCSV
    T = table(D.id,D.x_um,D.y_um,D.BORMatched,D.SF_Basal, ...
        D.SF_PrismaticA,D.SF_PyramidalCA,D.SF_MaxBasalOrPrismatic, ...
        D.CAxisToLoad_deg,D.BasalNormalStressOverSigma, ...
        'VariableNames',{'EBSD_ID','X_um','Y_um','BOR_matched', ...
        'SF_basal','SF_prismaticA','SF_pyramidalCA','SF_max_a', ...
        'CAxisToLoad_deg','BasalNormalStressOverSigma'});
    writetable(T,fullfile(folder,'17_Schmid_Alpha_PerPixel.csv'));
end
fprintf('Schmid analysis saved in: %s\n',folder);
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
    title(sprintf('%s | load axis [%g %g %g]',ttl,xyz),'Interpreter','none');
    % Directional arrow in figure coordinates, horizontally left for [-1,0,0].
    % Annotation orientation is illustrative, not a stress-vector simulation.
    d = xyz(1:2);
    if norm(d)>1e-10
        try
            d = d/norm(d);
            annotation('textarrow',[0.30 0.30+0.105*d(1)], ...
                [0.87 0.87+0.09*d(2)],'String','Load', ...
                'Color','black','LineWidth',1.4,'FontSize',10);
        catch MEarrow
            warning('Load-arrow annotation failed: %s',MEarrow.message);
        end
    end
    saveFigure(folder,filename);
catch ME
    warning('Could not save %s: %s',filename,ME.message);
end
end
