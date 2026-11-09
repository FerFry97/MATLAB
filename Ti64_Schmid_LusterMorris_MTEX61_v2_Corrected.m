%% Ti64_Schmid_LusterMorris_MTEX61_v2_Corrected.m
% Ti-6Al-4V / alpha-HCP: Schmid factor, Luster-Morris m', c-axis and
% orientation-based soft/hard CANDIDATES. MTEX 6.x / MATLAB R2019b+.
%
% v2 fixes:
%  * Excludes NOT-INDEXED pseudo-grains produced by calcGrains/minPixel;
%    excludes any indexed grain without finite orientation-based outputs.
%  * Computes one value per unique alpha/alpha grain pair, but weights
%    frequency distributions by total 2D boundary TRACE LENGTH (um).
%  * Adds c-axis/load-angle map, orientation-based soft/hard candidate map,
%    and candidate soft-hard interface map.
%  * Optional USER-SUPPLIED MTR polygon or interactive ROI. Within/outside
%    MTR comparisons are NOT generated without an explicitly defined ROI.
%  * Exports grain / boundary-pair tables and a reproducible analysis log.
%
% SCIENTIFIC LIMITS:
% Schmid and m' are geometric orientation descriptors (not measured slip).
% 'Soft' and 'Hard' below mean ORIENTATION-BASED CANDIDATES ONLY; the
% classification is NOT a strength measurement or a simulation.
% m' does NOT calculate local stress, load shedding, dwell-fatigue lifetime,
% or transmission events. EBSD 2D does not measure boundary plane normals.
% Confirm the CTF reference frame with your previous parent-beta script.
%
% MTEX docs:
% https://mtex-toolbox.github.io/GrainReconstruction.html
% https://mtex-toolbox.github.io/SlipTransmission.html
% https://mtex-toolbox.github.io/grainBoundary.segLength.html

%% ======================= USER OPTIONS ================================
opt.filename                 = '';          % select .ctf when empty
opt.phaseName                = '';          % detect HCP alpha or ask
opt.loadXYZ                  = [-1 0 0];    % sample-frame loading axis
opt.segAngleDeg              = 5;           % reconstruct alpha grains
opt.minPixels                = 5;           % minimum valid indexed grain
opt.smoothSteps              = 0;           % 0: RAW boundary lengths (recommended)
opt.showFigures              = true;
opt.boundaryLineWidth        = 1.8;
opt.outputDir                = '';          % alongside the input CTF

% PROXY classification for orientations, NOT actual local mechanical state
opt.softSchmidMin            = 0.42;  % max(m_basal, m_prismatic) >= this
opt.hardCAxisMaxDeg          = 15;    % c-axis roughly parallel to load
opt.hardSchmidMax            = 0.20;  % max(m_basal, m_prismatic) <= this

% MTR manual area: N-by-2 [sample X, sample Y] polygon (in CTF units).
% Option A: set polygon coordinates here; option B: set interactiveROI=true
% and click >=3 polygon vertices, then press Enter. Just Enter skips.
% NO automatic MTR identification or unverifiable MTR label is performed.
opt.mtrPolygonXY             = [];
opt.interactiveROI           = true;
%% ====================================================================

if exist('EBSD','class') ~= 8 || exist('slipSystem','class') ~= 8
    error('MTEX not found on the MATLAB path. Run startup_mtex first.');
end
if ~isnumeric(opt.loadXYZ) || numel(opt.loadXYZ)~=3 || ...
        any(~isfinite(opt.loadXYZ)) || norm(opt.loadXYZ) <= eps
    error('opt.loadXYZ must be a finite, nonzero three-component vector.');
end
if opt.minPixels<1 || opt.minPixels~=round(opt.minPixels)
    error('opt.minPixels must be a positive integer.');
end
if opt.softSchmidMin<0 || opt.softSchmidMin>0.5 || ...
        opt.hardSchmidMax<0 || opt.hardSchmidMax>0.5 || ...
        opt.hardCAxisMaxDeg<0 || opt.hardCAxisMaxDeg>90
    error('Orientation-based soft/hard thresholds are invalid.');
end

if isempty(opt.filename)
    [shortName, sourceDir] = uigetfile({'*.ctf','Oxford EBSD (*.ctf)'}, ...
        'Select an Oxford CTF map of Ti-6Al-4V');
    if isequal(shortName,0), error('Selection cancelled.'); end
    fileName = fullfile(sourceDir,shortName);
else
    fileName = opt.filename;
    [sourceDir,~,~] = fileparts(fileName);
    if isempty(sourceDir), sourceDir=pwd; end
end
if exist(fileName,'file')~=2, error('Input file not found: %s',fileName); end
if isempty(opt.outputDir)
    outDir=fullfile(sourceDir,'Ti64_Schmid_LusterMorris_v2_Output');
else
    outDir=opt.outputDir;
end
if exist(outDir,'dir')~=7, mkdir(outDir); end

L = vector3d(opt.loadXYZ(1),opt.loadXYZ(2),opt.loadXYZ(3));
L = L ./ norm(L);
vis='on'; if ~opt.showFigures, vis='off'; end

fprintf('\n======= Ti64 SCHMID / LUSTER-MORRIS v2 =======\n');
fprintf('CTF: %s\nLoad axis: [%g %g %g] (sample coordinates)\n', ...
    fileName,opt.loadXYZ);
fprintf('Grain threshold = %g deg; minPixels = %d\n', ...
    opt.segAngleDeg,opt.minPixels);

%% 1) Import only alpha indexed data
% Match this option to the reference frame of your other MTEX programs.
ebAll=EBSD.load(fileName,'interface','ctf', ...
    'convertEuler2SpatialReferenceFrame');
allNames=ebAll.mineralList;
if ischar(allNames), allNames=cellstr(allNames); end
if isstring(allNames), allNames=cellstr(allNames); end
keep=~strcmpi(allNames,'notIndexed') & ~cellfun(@isempty,allNames);
phaseNames=allNames(keep);
if isempty(phaseNames), error('No indexed crystalline phases in CTF.'); end

if isempty(opt.phaseName)
    s=lower(string(phaseNames));
    iAlpha=contains(s,'alpha') | contains(s,'hcp') | contains(s,char(945));
    if nnz(iAlpha)==1
        phaseName=phaseNames{find(iAlpha,1)};
    elseif numel(phaseNames)==1
        phaseName=phaseNames{1};
        warning('One phase found. VERIFY that %s is HCP alpha.',phaseName);
    else
        fprintf('Indexed phases:\n');
        for k=1:numel(phaseNames), fprintf('%d) %s\n',k,phaseNames{k}); end
        k=input('Choose HCP alpha phase NUMBER: ');
        if isempty(k)||~isscalar(k)||~isfinite(k)||k<1|| ...
                k>numel(phaseNames)||k~=round(k)
            error('Invalid alpha selection. Configure opt.phaseName.');
        end
        phaseName=phaseNames{k};
    end
else
    phaseName=opt.phaseName;
end

ebAlpha=ebAll(phaseName);
if isempty(ebAlpha), error('Alpha phase has no indexed data.'); end
cs=ebAlpha.CS;
fprintf('Selected phase: %s\n',phaseName);
disp(cs);

%% 2) Reconstruct; EXCLUDE notIndexed pseudo-grains, not only small pixels
[grainsAll,ebAlpha]=calcGrains(ebAlpha,'angle', ...
    opt.segAngleDeg*degree,'minPixel',opt.minPixels);
numAll=length(grainsAll);
% IMPORTANT: calcGrains(minPixel) can return notIndexed grains with 1-2
% pixels and NaN meanOrientation. The old script accidentally exported them.
grains=grainsAll('indexed');
fprintf('Grain objects: all=%d, indexed=%d; pseudo-grains excluded=%d\n', ...
    numAll,length(grains),numAll-length(grains));
if isempty(grains), error('No indexed alpha grains remain.'); end

% No smoothing is applied to grain boundaries by default. If enabled,
% weighted lengths refer to the smoothed *2D trace*, not the original trace.
if opt.smoothSteps>0
    warning(['Boundary smoothing changes trace lengths. For reproducible ' ...
             'raw-length weighting, opt.smoothSteps=0 is recommended.']);
    grains=smoothBoundary(grains,opt.smoothSteps);
end

%% 3) Grain geometry and Schmid (valid grains only)
oriG=grains.meanOrientation;
sBas=slipSystem.basal(cs).symmetrise('antipodal');
sPri=slipSystem.prismaticA(cs).symmetrise('antipodal');
sPyr=slipSystem.pyramidalCA(cs).symmetrise('antipodal');

[SFbas,iBas]=max(abs(sBas.SchmidFactor(inv(oriG)*L)),[],2);
[SFpri,iPri]=max(abs(sPri.SchmidFactor(inv(oriG)*L)),[],2);
[SFpyr,~]=max(abs(sPyr.SchmidFactor(inv(oriG)*L)),[],2);
cAxis=Miller(0,0,0,1,cs,'hkl');
thetaC=angle(oriG*cAxis,L)./degree;
thetaC=min(thetaC(:),180-thetaC(:));

finiteGrain=isfinite(SFbas(:)) & isfinite(SFpri(:)) & ...
    isfinite(SFpyr(:)) & isfinite(thetaC(:)) & ...
    grains.numPixel(:)>=opt.minPixels;
if any(~finiteGrain)
    warning('Dropping %d indexed grains with invalid orientations/size.', ...
        nnz(~finiteGrain));
    grains=grains(finiteGrain);
    oriG=grains.meanOrientation;
    [SFbas,iBas]=max(abs(sBas.SchmidFactor(inv(oriG)*L)),[],2);
    [SFpri,iPri]=max(abs(sPri.SchmidFactor(inv(oriG)*L)),[],2);
    [SFpyr,~]=max(abs(sPyr.SchmidFactor(inv(oriG)*L)),[],2);
    thetaC=angle(oriG*cAxis,L)./degree;
    thetaC=min(thetaC(:),180-thetaC(:));
end
if isempty(grains), error('No valid grains remain.'); end
assert(all(isfinite([SFbas(:);SFpri(:);SFpyr(:);thetaC(:)])), ...
    'Invalid values remain in Schmid/c-axis arrays.');

% These labels are GEOMETRIC PROXIES, NOT measured strength or creep rates.
maxEasySF=max([SFbas(:),SFpri(:)],[],2);
isHard=thetaC<=opt.hardCAxisMaxDeg & maxEasySF<=opt.hardSchmidMax;
isSoft=maxEasySF>=opt.softSchmidMin & ~isHard;
classNum=zeros(length(grains),1); % 0=Other, 1=Soft candidate, 2=Hard candidate
classNum(isSoft)=1;
classNum(isHard)=2;
className=repmat("Other",length(grains),1);
className(isSoft)="SoftCandidate";
className(isHard)="HardCandidate";

%% 4) Optional MTR ROI in sample coordinates
roiPoly=opt.mtrPolygonXY;
if ~isempty(roiPoly)
    if ~isnumeric(roiPoly) || size(roiPoly,2)~=2 || ...
            size(roiPoly,1)<3 || any(~isfinite(roiPoly(:)))
        error('opt.mtrPolygonXY must be finite N-by-2 with >=3 vertices.');
    end
end
if isempty(roiPoly) && opt.interactiveROI
    if ~opt.showFigures
        warning('interactiveROI requires visible figures; skipping ROI.');
    else
        hPick=figure('Color','w','Name','Choose MTR ROI','Visible','on');
        plot(grains,grains.meanOrientation,'ipfDirection',L,'figSize','large');
        title('IPF toward load axis: click MTR outline, then press ENTER');
        fprintf(['MTR ROI: CLICK >= 3 vertices around the identified MTR, ' ...
                 'then press ENTER. Enter alone skips ROI.\n']);
        [xp,yp]=ginput;
        if numel(xp)>=3, roiPoly=[xp(:),yp(:)]; end
        close(hPick);
    end
end
cg=grains.centroid;
roiDefined=~isempty(roiPoly);
isInMTR=false(length(grains),1);
regionName=repmat("NotDefined",length(grains),1);
if roiDefined
    isInMTR=inpolygon(cg.x(:),cg.y(:),roiPoly(:,1),roiPoly(:,2));
    regionName(:)="OutsideMTR";
    regionName(isInMTR)="InsideMTR";
    fprintf(['WARNING: MTR classification is based on grain CENTROIDS in ' ...
             'a 2D manually defined ROI. Inspect crossing grains.\n']);
else
    fprintf('No MTR polygon: MTR-versus-outside comparison disabled.\n');
end

%% 5) All valid internal alpha/alpha boundaries; make unique grain pairs
allGB=grains.boundary('indexed');
if isempty(allGB), error('No indexed alpha/alpha boundary segments.'); end
segIDs=allGB.grainId;
keepSeg=ismember(segIDs(:,1),grains.id) & ...
    ismember(segIDs(:,2),grains.id) & ...
    segIDs(:,1)>0 & segIDs(:,2)>0 & ...
    segIDs(:,1)~=segIDs(:,2);
allGB=allGB(keepSeg);
if isempty(allGB), error('No two-sided alpha/alpha boundaries remain.'); end
segLen=double(allGB.segLength(:));
validSeg=isfinite(segLen)&segLen>0;
if any(~validSeg)
    warning('Removing %d nonpositive/invalid boundary-length segments.', ...
        nnz(~validSeg));
    allGB=allGB(validSeg);
    segLen=segLen(validSeg);
end
if isempty(allGB), error('No valid alpha/alpha boundary trace remains.'); end

gB=allGB;
[pairIDs,~,pairIdx]=unique(sort(gB.grainId,2),'rows');
nPairs=size(pairIDs,1);
ind1=grains.id2ind(pairIDs(:,1));
ind2=grains.id2ind(pairIDs(:,2));
if any(~isfinite(ind1))||any(~isfinite(ind2))|| ...
        any(ind1<1)||any(ind2<1)
    error('A grain-boundary pair references missing indexed grains.');
end
pairLen=accumarray(pairIdx,segLen,[nPairs,1],@sum,0);
pairSegments=accumarray(pairIdx,1,[nPairs,1],@sum,0);
fprintf('VALID alpha grains: %d; segments: %d; pairs: %d\n', ...
    length(grains),length(gB),nPairs);

%% 6) Schmid-selected and maximum geometric Luster-Morris m'
sBasChosen=oriG.*sBas(iBas);
sPriChosen=oriG.*sPri(iPri);
mSelBas=abs(mPrime(sBasChosen(ind1),sBasChosen(ind2)));
mSelPri=abs(mPrime(sPriChosen(ind1),sPriChosen(ind2)));
mBestBas=bestPairwiseMPrime(oriG(ind1),oriG(ind2),sBas);
mBestPri=bestPairwiseMPrime(oriG(ind1),oriG(ind2),sPri);
mSelBas=min(1,max(0,mSelBas(:)));
mSelPri=min(1,max(0,mSelPri(:)));
mBestBas=min(1,max(0,mBestBas(:)));
mBestPri=min(1,max(0,mBestPri(:)));
misDeg=angle(oriG(ind1),oriG(ind2))./degree;
misDeg=misDeg(:);

if any(~isfinite([mSelBas;mSelPri;mBestBas;mBestPri;misDeg;pairLen]))
    error('Nonfinite boundary quantity detected; check orientations.');
end
if any(mBestBas+1e-9<mSelBas) || any(mBestPri+1e-9<mSelPri)
    warning(['Geometric maximum is lower than load-selected m'' for one or ' ...
             'more pairs. Check slip-system symmetry and candidate enumeration.']);
end

% Pair classification: 1=soft-hard, 2=soft-soft, 3=hard-hard, 0=other
isSoftHard=(isSoft(ind1)&isHard(ind2)) | (isHard(ind1)&isSoft(ind2));
classPair=repmat("OtherPair",nPairs,1);
classPair(isSoftHard)="SoftHardCandidate";
classPair(isSoft(ind1)&isSoft(ind2))="SoftSoftCandidate";
classPair(isHard(ind1)&isHard(ind2))="HardHardCandidate";

regionPair=repmat("NotAssessed",nPairs,1);
if roiDefined
    bothIn=isInMTR(ind1)&isInMTR(ind2);
    bothOut=~isInMTR(ind1)&~isInMTR(ind2);
    regionPair(:)="MTR_Interface";
    regionPair(bothIn)="MTR_Inside";
    regionPair(bothOut)="MTR_Outside";
end

%% 7) Figures (01-09: preserve old analysis; 10-12: new)
h=plotGrainQuantity(grains,SFbas,'Schmid |m|: basal <a>', ...
    'Basal |Schmid|',[0 .5],vis);
saveFigure(h,outDir,'01_Schmid_Basal');
h=plotGrainQuantity(grains,SFpri,'Schmid |m|: prismatic <a>', ...
    'Prismatic |Schmid|',[0 .5],vis);
saveFigure(h,outDir,'02_Schmid_Prismatic');
h=plotGrainQuantity(grains,SFpyr,'Schmid |m|: pyramidal <c+a>', ...
    'Pyramidal <c+a> |Schmid|',[0 .5],vis);
saveFigure(h,outDir,'03_Schmid_PyramidalCA');

h=plotBoundaryQuantity(grains,gB,mSelBas(pairIdx), ...
    'm'' basal: Schmid-selected','Selected basal m''',[0 1], ...
    opt.boundaryLineWidth,vis);
saveFigure(h,outDir,'04_mPrime_Selected_Basal');
h=plotBoundaryQuantity(grains,gB,mSelPri(pairIdx), ...
    'm'' prismatic: Schmid-selected','Selected prismatic m''',[0 1], ...
    opt.boundaryLineWidth,vis);
saveFigure(h,outDir,'05_mPrime_Selected_Prismatic');
h=plotBoundaryQuantity(grains,gB,mBestBas(pairIdx), ...
    'm'' basal: geometric maximum','Maximum basal m''',[0 1], ...
    opt.boundaryLineWidth,vis);
saveFigure(h,outDir,'06_mPrime_Best_Basal');
h=plotBoundaryQuantity(grains,gB,mBestPri(pairIdx), ...
    'm'' prismatic: geometric maximum','Maximum prismatic m''',[0 1], ...
    opt.boundaryLineWidth,vis);
saveFigure(h,outDir,'07_mPrime_Best_Prismatic');

edges=0:0.05:1;
centers=(edges(1:end-1)+edges(2:end))/2;
h=figure('Color','w','Name','Length-weighted m'' distribution','Visible',vis);
y1=weightedBins(mSelBas,pairLen,edges);
y2=weightedBins(mBestBas,pairLen,edges);
plot(centers,y1,'LineWidth',1.8); hold on;
plot(centers,y2,'LineWidth',1.8); hold off;
xlabel('Luster-Morris m''');
ylabel('Fraction of total alpha/alpha boundary TRACE LENGTH');
legend({'Basal Schmid-selected','Basal geometric maximum'},'Location','best');
title('Basal m'' distributions, weighted by 2D trace length');
xlim([0 1]); grid on;
saveFigure(h,outDir,'08_mPrime_LengthWeighted_Histogram');

h=figure('Color','w','Name','m'' versus misorientation','Visible',vis);
scatter(misDeg,mBestBas,max(10,10*sqrt(pairLen./median(pairLen))), ...
    mSelBas,'filled');
xlabel('alpha/alpha grain-mean misorientation (degrees)');
ylabel('Geometric maximum basal m''');
title('Each dot = one grain pair; marker area scales with boundary length');
colormap(parula(256)); caxis([0 1]);
cb=colorbar; ylabel(cb,'Schmid-selected basal m''');
ylim([0 1]); grid on;
saveFigure(h,outDir,'09_mPrime_vs_Misorientation');

h=plotGrainQuantity(grains,thetaC, ...
    'c-axis declination angle vs loading axis', ...
    'Acute angle [0001] to load (deg)',[0 90],vis);
saveFigure(h,outDir,'10_CAxis_Load_Angle');

h=figure('Color','w','Name','Orientation-based soft/hard candidates', ...
    'Visible',vis);
plot(grains,classNum,'figSize','large'); hold on;
plot(grains.boundary,'lineWidth',0.3,'lineColor',[0.2 0.2 0.2]);
hold off;
colormap([0.73 0.73 0.73; 0.20 0.62 0.82; 0.91 0.35 0.24]);
caxis([-0.5 2.5]); cb=colorbar; cb.Ticks=[0 1 2];
cb.TickLabels={'Other','Soft candidate','Hard candidate'};
title('Orientation-based candidates (NOT actual stress/strain)');
axis equal tight;
saveFigure(h,outDir,'11_SoftHard_Orientation_Candidates');

h=figure('Color','w','Name','Candidate soft-hard alpha/alpha interfaces', ...
    'Visible',vis);
plot(grains,'FaceColor',[.89 .89 .89],'figSize','large'); hold on;
plot(gB,'lineColor',[.73 .73 .73],'lineWidth',0.5);
if any(isSoftHard)
    plot(gB(isSoftHard(pairIdx)),'lineColor',[.82 .1 .12], ...
        'lineWidth',max(2,opt.boundaryLineWidth));
end
hold off;
axis equal tight;
title('Soft/hard candidate interfaces (red); proxy classification only');
saveFigure(h,outDir,'12_SoftHard_Boundary_Candidates');

%% 8) MTR comparison ONLY with an explicitly defined user ROI
if roiDefined
    h=figure('Color','w','Name','MTR vs outside: weighted basal m''', ...
        'Visible',vis);
    maskIn=(regionPair=="MTR_Inside");
    maskOut=(regionPair=="MTR_Outside");
    yi=weightedBins(mBestBas(maskIn),pairLen(maskIn),edges);
    yo=weightedBins(mBestBas(maskOut),pairLen(maskOut),edges);
    plot(centers,yi,'LineWidth',1.8); hold on;
    plot(centers,yo,'LineWidth',1.8); hold off;
    legend({'Both grain centroids inside MTR', ...
            'Both grain centroids outside MTR'},'Location','best');
    xlabel('Maximum geometric basal m''');
    ylabel('Fraction of trace length IN EACH group');
    xlim([0 1]); grid on;
    title('Manual-ROI MTR comparison (interface pairs excluded)');
    saveFigure(h,outDir,'13_MTR_vs_Outside_Basal_mPrime');

    h=figure('Color','w','Name','Manual MTR polygon','Visible',vis);
    plot(grains,thetaC,'figSize','large'); hold on;
    xx=[roiPoly(:,1);roiPoly(1,1)];
    yy=[roiPoly(:,2);roiPoly(1,2)];
    plot(xx,yy,'k-','LineWidth',2); hold off;
    colormap(parula(256)); caxis([0 90]);
    mtexColorbar('title','c-axis angle (deg)');
    title('Manual ROI overlay (membership by grain centroid)');
    saveFigure(h,outDir,'14_MTR_Manual_ROI');
end

%% 9) CSV exports: valid grains only + one row per UNIQUE alpha/alpha pair
cg=grains.centroid;
grainTable=table(grains.id(:),cg.x(:),cg.y(:),grains.numPixel(:), ...
    thetaC(:),SFbas(:),SFpri(:),SFpyr(:),maxEasySF(:), ...
    className,regionName, ...
    'VariableNames',{'GrainID','CentroidX','CentroidY','NumPixels', ...
    'Caxis_Load_Angle_deg','Schmid_Basal','Schmid_Prismatic', ...
    'Schmid_PyramidalCA','Max_Schmid_BasalOrPrismatic', ...
    'OrientationClass_Proxy','MTR_GrainCentroid_Region'});
writetable(grainTable,fullfile(outDir,'Grain_Schmid_Corrected.csv'));

pairTable=table(pairIDs(:,1),pairIDs(:,2),pairLen(:),pairSegments(:), ...
    misDeg,mSelBas,mSelPri,mBestBas,mBestPri, ...
    thetaC(ind1),thetaC(ind2),classPair,regionPair, ...
    'VariableNames',{'GrainID_1','GrainID_2','BoundaryTraceLength_um', ...
    'NumBoundarySegments','Misorientation_deg','mPrime_Selected_Basal', ...
    'mPrime_Selected_Prismatic','mPrime_Best_Basal', ...
    'mPrime_Best_Prismatic','Caxis_Angle_Grain1_deg', ...
    'Caxis_Angle_Grain2_deg','SoftHardPair_Proxy', ...
    'MTR_BoundaryRegion_ManualROI'});
writetable(pairTable,fullfile(outDir,'Boundary_mPrime_Corrected.csv'));

if roiDefined
    labels=["MTR_Inside";"MTR_Outside";"MTR_Interface"];
    nPairsG=zeros(3,1); totalLength=zeros(3,1);
    meanBestBas=nan(3,1); medBestBas=nan(3,1);
    meanBestPri=nan(3,1);
    fracSoftHardLength=nan(3,1);
    for k=1:3
        ix=regionPair==labels(k);
        nPairsG(k)=nnz(ix);
        totalLength(k)=sum(pairLen(ix));
        if any(ix)
            meanBestBas(k)=sum(mBestBas(ix).*pairLen(ix))/totalLength(k);
            medBestBas(k)=weightedMedian(mBestBas(ix),pairLen(ix));
            meanBestPri(k)=sum(mBestPri(ix).*pairLen(ix))/totalLength(k);
            fracSoftHardLength(k)=sum(pairLen(ix & isSoftHard))/totalLength(k);
        end
    end
    roiSummary=table(labels,nPairsG,totalLength,meanBestBas,medBestBas, ...
        meanBestPri,fracSoftHardLength,'VariableNames', ...
        {'Region','UniqueGrainPairs','TotalTraceLength_um', ...
        'LengthWeightedMean_BestBasal_mPrime', ...
        'LengthWeightedMedian_BestBasal_mPrime', ...
        'LengthWeightedMean_BestPrismatic_mPrime', ...
        'FractionLength_SoftHardCandidate'});
    writetable(roiSummary,fullfile(outDir,'MTR_ROI_Comparison.csv'));
    writematrix(roiPoly,fullfile(outDir,'MTR_ROI_vertices_XY.csv'));
end

%% 10) Numerical summary and parameter log
fid=fopen(fullfile(outDir,'Analysis_Parameters_And_Statistics.txt'),'w');
if fid<0, warning('Could not write Analysis_Parameters_And_Statistics.txt'); end
reportLine(fid,'Input file: %s\n',fileName);
reportLine(fid,'Alpha phase: %s\n',phaseName);
reportLine(fid,'Import: convertEuler2SpatialReferenceFrame\n');
reportLine(fid,'Load vector (not necessarily verified physical load): [%g %g %g]\n',opt.loadXYZ);
reportLine(fid,'Grain threshold (deg): %.2f; minPixels: %d; smoothing: %g\n', ...
    opt.segAngleDeg,opt.minPixels,opt.smoothSteps);
reportLine(fid,'Valid indexed alpha grains: %d; pseudo-grains excluded: %d\n', ...
    length(grains),numAll-length(grains));
reportLine(fid,'Valid alpha/alpha boundary segments: %d; unique pairs: %d\n', ...
    length(gB),nPairs);
reportLine(fid,'Total 2D alpha/alpha boundary trace length (um): %.5f\n',sum(pairLen));
reportLine(fid,'Schmid (valid alpha grains): basal mean=%.4f; prism mean=%.4f; pyr<c+a> mean=%.4f\n', ...
    mean(SFbas),mean(SFpri),mean(SFpyr));
reportLine(fid,'m'' basal selected: UNWEIGHTED pair median=%.4f, TRACE-WEIGHTED median=%.4f\n', ...
    median(mSelBas),weightedMedian(mSelBas,pairLen));
reportLine(fid,'m'' basal max: UNWEIGHTED pair median=%.4f, TRACE-WEIGHTED median=%.4f\n', ...
    median(mBestBas),weightedMedian(mBestBas,pairLen));
reportLine(fid,'m'' prism selected: UNWEIGHTED median=%.4f, TRACE-WEIGHTED median=%.4f\n', ...
    median(mSelPri),weightedMedian(mSelPri,pairLen));
reportLine(fid,'m'' prism max: UNWEIGHTED median=%.4f, TRACE-WEIGHTED median=%.4f\n', ...
    median(mBestPri),weightedMedian(mBestPri,pairLen));
reportLine(fid,'Soft proxy >= %.3f (max basal/prismatic Schmid).\n',opt.softSchmidMin);
reportLine(fid,'Hard proxy: thetaC <= %.2f deg AND max basal/prism Schmid <= %.3f.\n', ...
    opt.hardCAxisMaxDeg,opt.hardSchmidMax);
reportLine(fid,'Soft candidate grains: %d; hard candidate grains: %d\n', ...
    nnz(isSoft),nnz(isHard));
reportLine(fid,'Candidate soft-hard boundary pairs: %d of %d (%.2f %% trace length).\n', ...
    nnz(isSoftHard),nPairs,100*sum(pairLen(isSoftHard))/sum(pairLen));
reportLine(fid,'MTR manual ROI defined: %d. If yes, grain centroids determine region.\n',roiDefined);
reportLine(fid,'Interpretation: NO local stress field, load shedding or measured slip is modeled.\n');
if fid>=0, fclose(fid); end
fprintf('\nDONE: %s\n',outDir);
fprintf('Indexed alpha grains=%d, unique alpha/alpha pairs=%d\n', ...
    length(grains),nPairs);
fprintf('m'' basal best, weighted median = %.3f\n', ...
    weightedMedian(mBestBas,pairLen));
fprintf('NOTE: MTR comparisons only if a manual polygon was supplied.\n');

%% ========================== LOCAL FUNCTIONS ============================
function result=bestPairwiseMPrime(ori1,ori2,sS)
    % All pair combinations, one maximum per grain-pair. Chunked by slip
    % variants (no nPairs x nPairs matrix is generated).
    n=length(ori1);
    result=zeros(n,1);
    for ia=1:length(sS)
        sA=ori1.*sS(ia);
        for ib=1:length(sS)
            sB=ori2.*sS(ib);
            candidate=abs(mPrime(sA,sB));
            result=max(result,candidate(:));
        end
    end
end

function mass=weightedBins(values,weights,edges)
    n=numel(edges)-1;
    mass=zeros(n,1);
    if isempty(values), return; end
    [~,~,id]=histcounts(values(:),edges);
    valid=id>0 & isfinite(values(:)) & ...
        isfinite(weights(:)) & weights(:)>0;
    if any(valid)
        mass=accumarray(id(valid),weights(valid),[n,1],@sum,0);
        mass=mass/sum(mass);
    end
end

function med=weightedMedian(values,weights)
    valid=isfinite(values(:))&isfinite(weights(:))&weights(:)>0;
    med=NaN;
    if ~any(valid), return; end
    v=values(valid); w=weights(valid);
    [v,order]=sort(v(:)); w=w(order);
    j=find(cumsum(w)>=sum(w)/2,1,'first');
    med=v(j);
end

function h=plotGrainQuantity(grains,values,name,barTitle,rng,vis)
    h=figure('Color','w','Name',name,'Visible',vis);
    plot(grains,values(:),'figSize','large'); hold on;
    plot(grains.boundary,'lineWidth',0.3,'lineColor',[.22 .22 .22]);
    hold off; colormap(parula(256)); caxis(rng);
    mtexColorbar('title',barTitle);
    title(name,'Interpreter','none'); axis equal tight;
end

function h=plotBoundaryQuantity(grains,gB,values,name,barTitle,rng,lw,vis)
    h=figure('Color','w','Name',name,'Visible',vis);
    plot(grains,'FaceColor',[.90 .90 .90],'figSize','large'); hold on;
    plot(gB,values(:),'lineWidth',lw); hold off;
    colormap(parula(256)); caxis(rng);
    mtexColorbar('title',barTitle);
    title(name,'Interpreter','none'); axis equal tight;
end

function saveFigure(h,folder,base)
    drawnow;
    savefig(h,fullfile(folder,[base '.fig']));
    try
        exportgraphics(h,fullfile(folder,[base '.png']),'Resolution',250);
    catch ME
        warning('exportgraphics failed for %s: %s; trying print',base,ME.message);
        print(h,fullfile(folder,[base '.png']),'-dpng','-r250');
    end
end

function reportLine(fid,varargin)
    fprintf(varargin{:});
    if fid>=0, fprintf(fid,varargin{:}); end
end
