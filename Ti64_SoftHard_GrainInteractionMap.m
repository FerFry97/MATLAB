function Ti64_SoftHard_GrainInteractionMap()
% Standalone soft-hard alpha grain map. Requires MTEX + original .ctf EBSD.
% Does NOT re-run parent-beta reconstruction and saves into a SEPARATE folder.
% In the MATLAB Command Window run: Ti64_SoftHard_GrainInteractionMap
% Authoring note: this code has not been runtime-tested in MATLAB/MTEX.
clc;
if exist('EBSD','class')~=8 || exist('slipSystem','class')~=8
    candidate=which('startup_mtex');
    if isempty(candidate)
        root=uigetdir(pwd,'Select the FULL MTEX installation folder');
        if isequal(root,0),error('MTEX startup cancelled.');end
    else
        root=fileparts(candidate);
    end
    if exist(fullfile(root,'tools','option_tools','setMTEXpref.m'),'file')~=2
        root=uigetdir(pwd,['The startup_mtex.m on your MATLAB path is not a full ' ...
            'installation. Select the actual MTEX root folder']);
        if isequal(root,0),error('A complete MTEX installation is required.');end
    end
    if exist(fullfile(root,'tools','option_tools','setMTEXpref.m'),'file')~=2
        error('The selected folder is not a complete MTEX installation.');
    end
    cdBefore=pwd;
    addpath(root);
    try
        cd(root); startup_mtex('noMenu');
    catch ME
        cd(cdBefore);rethrow(ME);
    end
    cd(cdBefore);
end
alphaPhase='Ti-Hex';
loadAxisXYZ=[-1 0 0];
nominalStressMPa=855; % purely assumed UNIFORM tension, not a FEM model
[fn,fp]=uigetfile({'*.ctf','Oxford AZtec CTF (*.ctf)'}, ...
    'Select Sample7.ctf or SAMPLE8_crop.ctf');
if isequal(fn,0),disp('No CTF chosen.');return;end
input=fullfile(fp,fn);
[~,nm]=fileparts(fn);
folder=fullfile(fp,[nm '_SoftHard_AlphaGrainMap']);
if ~isfolder(folder),mkdir(folder);end
fprintf('Loading measured EBSD alpha from %s\n',input);
ebsd=EBSD.load(input,'interface','ctf','convertEuler2SpatialReferenceFrame');
alpha=ebsd(alphaPhase);
if isempty(alpha),error('No indexed alpha phase named %s in this CTF.',alphaPhase);end
opt=struct('alphaPhase',alphaPhase,'loadAxisXYZ',loadAxisXYZ,...
 'nominalStressMPa',nominalStressMPa,'grainAngleDeg',5,...
 'minGrainPixels',30,'softSF',0.40,'softMinCAngleDeg',50,...
 'hardSF',0.30,'hardMaxCAngleDeg',30);
OUT=softHardGrainMap(alpha,alpha.CS,opt,folder); %#ok<NASGU>
fprintf('\nDONE: maps 22-24 and two CSV files in:\n%s\n',folder);
end

function OUT = softHardGrainMap(alphaMeasured,CSalpha,opt,folder)
% SOFT/HARD GRAIN INTERACTION SCREENING (MEASURED alpha orientations ONLY).
% Does NOT depend on prior-beta reconstruction; does NOT solve load shedding.
% Grain units, not pixels: one mean orientation per segmented alpha grain.
% Axis is a uniaxial tensile loading axis in the specimen/EBSD frame.
% "soft" and "hard" are geometric screening categories, NOT proven mechanics.
% CRSS and neighbor constraints must be known for an actual plasticity model.

if ~isfolder(folder), mkdir(folder); end
nPix0 = numel(alphaMeasured.id);
xyz=double(opt.loadAxisXYZ(:)');
assert(numel(xyz)==3 && all(isfinite(xyz)) && norm(xyz)>0, ...
    'SHopt.loadAxisXYZ must be a valid nonzero vector of length 3.');
xyz=xyz/norm(xyz);
assert(opt.grainAngleDeg>0 && opt.grainAngleDeg<30,'Choose 0 < grainAngleDeg < 30');
assert(opt.minGrainPixels>=1,'minGrainPixels must be >=1');
assert(opt.softSF>opt.hardSF,'softSF must exceed hardSF');
loadDir=vector3d(xyz(1),xyz(2),xyz(3));

% MTEX 6.x/7.x supported form: returns grain2d and a copy with grainId.
% Avoid silently modifying the full EBSD/reconstruction data upstream.
[gs,alphaWithIds] = calcGrains(alphaMeasured,'angle',opt.grainAngleDeg*degree,...
    'minPixel',opt.minGrainPixels);
% When called with an alpha-only EBSD, gs can still contain unindexed areas;
% take alpha only, using the supplied phase name where possible.
try
    gs=gs(opt.alphaPhase);
catch
    error('Cannot filter alpha grains by phase %s. Check alphaPhase.',opt.alphaPhase);
end
if isempty(gs), error('No alpha grains after segmentation. Reduce minGrainPixels.'); end
ids=double(gs.id(:)); n=length(gs);
ori=gs.meanOrientation;
ln = inv(ori)*loadDir;
ssB=slipSystem.basal(CSalpha).symmetrise('antipodal');
ssP=slipSystem.prismaticA(CSalpha).symmetrise('antipodal');
ssCA=slipSystem.pyramidalCA(CSalpha).symmetrise('antipodal');
mB=max(abs(ssB.SchmidFactor(ln)),[],2); mB=mB(:);
mP=max(abs(ssP.SchmidFactor(ln)),[],2); mP=mP(:);
mCA=max(abs(ssCA.SchmidFactor(ln)),[],2); mCA=mCA(:);
mEasy=max(mB,mP);
cax=Miller(0,0,0,1,CSalpha,'UVTW');
cang=angle(ori*cax,loadDir)/degree;
cang=cang(:); cang=min(cang,180-cang);
sigN=cosd(cang).^2;
% Screening definitions: SHOPT PARAMETRIC, not universal material thresholds.
% Hard: c axis near load (high basal-plane normal traction) AND low easy-<a> SF.
% Soft: high easy-<a> SF AND c axis away from load axis.
classCode=2*ones(n,1);  % 1 hard, 2 intermediate, 3 soft
hard=(mEasy<=opt.hardSF) & (cang<=opt.hardMaxCAngleDeg);
soft=(mEasy>=opt.softSF) & (cang>=opt.softMinCAngleDeg);
classCode(hard)=1; classCode(soft)=3;
classNames={'Hard-oriented';'Intermediate';'Soft-oriented'};
numPixel=double(gs.grainSize(:));
areaUm2=double(area(gs));areaUm2=areaUm2(:);
ecdUm=2*sqrt(areaUm2/pi);
% Provisional orientation-dependent stress components under UNIFORM tension.
% Stress is not computed from rotation or FE simulations.
sigma=opt.nominalStressMPa;
if ~isfinite(sigma) || sigma<=0, sigma=NaN; end
Tgr=table(ids,numPixel,areaUm2,ecdUm,cang,mB,mP,mCA,mEasy,sigN, ...
    sigma*sigN,classCode,classNames(classCode), ...
    'VariableNames',{'GrainID','IndexedPixels','Area_um2','ECD_um',...
    'CAxisToLoad_deg','SF_Basal','SF_PrismaticA','SF_PyramidalCA',...
    'SF_MaxEasyA','BasalNormalStressOverSigma','BasalNormalStress_MPa',...
    'ClassCode','Class'});
writetable(Tgr,fullfile(folder,'22_SoftHard_AlphaGrains.csv'));

% Grain-boundary segment pairs (MTEX physical neighbor geometry).
% Only alpha-alpha boundaries between TWO indexed and retained grains.
B=gs.boundary;
gid=double(B.grainId);
[tf1,loc1]=ismember(gid(:,1),ids);
[tf2,loc2]=ismember(gid(:,2),ids);
valid=tf1 & tf2 & loc1>0 & loc2>0 & loc1~=loc2;
B=B(valid);
loc1=loc1(valid);loc2=loc2(valid);
segLen=double(B.segLength(:));
segIds=[ids(loc1),ids(loc2)];
segContrast=abs(mEasy(loc1)-mEasy(loc2));
segMis=angle(ori(loc1),ori(loc2))/degree;segMis=segMis(:);
segSH=(classCode(loc1)==1 & classCode(loc2)==3) | ...
      (classCode(loc1)==3 & classCode(loc2)==1);
if isempty(segLen)
    warning('No indexable alpha-alpha boundaries; check segmentation and indexing.');
    Tpair=table();
    nContact=0;contactLen=0;
else
    ordered=sort(segIds,2);
    [pairIDs,~,ic]=unique(ordered,'rows');
    nPairs=size(pairIDs,1);
    [~,iA]=ismember(pairIDs(:,1),ids);
    [~,iB]=ismember(pairIDs(:,2),ids);
    lengthUm=accumarray(ic,segLen,[nPairs,1],@sum,0);
    nSegments=accumarray(ic,1,[nPairs,1],@sum,0);
    dM=abs(mEasy(iA)-mEasy(iB));
    dTheta=abs(cang(iA)-cang(iB));
    misAng=angle(ori(iA),ori(iB))/degree;misAng=misAng(:);
    SH=(classCode(iA)==1 & classCode(iB)==3) | ...
       (classCode(iA)==3 & classCode(iB)==1);
    % purely geometric screening index, NOT a measured stress concentration
    proxyScore=dM .* max(sigN(iA),sigN(iB));
    Tpair=table(pairIDs(:,1),pairIDs(:,2),nSegments,lengthUm, ...
       classNames(classCode(iA)),classNames(classCode(iB)), ...
       mEasy(iA),mEasy(iB),cang(iA),cang(iB), ...
       mB(iA),mB(iB),mP(iA),mP(iB), ...
       dM,dTheta,misAng,sigN(iA),sigN(iB),SH,proxyScore, ...
       'VariableNames',{'GrainA','GrainB','BoundarySegments','ContactLength_um',...
       'ClassA','ClassB','SF_EasyA','SF_EasyB',...
       'CAxis_deg_A','CAxis_deg_B','SF_Basal_A','SF_Basal_B',...
       'SF_Prismatic_A','SF_Prismatic_B','DeltaSF_EasyA',...
       'DeltaCAxis_deg','Misorientation_deg','BasalNormalRatio_A',...
       'BasalNormalRatio_B','SoftHardContact','GeometricScreeningIndex'});
    Tpair=sortrows(Tpair,'GeometricScreeningIndex','descend');
    nContact=nnz(SH); contactLen=sum(lengthUm(SH));
end
writetable(Tpair,fullfile(folder,'23_SoftHard_NeighborPairs.csv'));

fprintf('\n--- SOFT-HARD GRAIN INTERACTION SCREEN ---\n');
fprintf('Alpha measured pixels: %d | segmented alpha grains: %d\n',nPix0,n);
fprintf('Grain angle %.1f deg; min grain pixels %d\n',opt.grainAngleDeg,opt.minGrainPixels);
fprintf('Hard-oriented grains: %d | soft-oriented grains: %d | intermediate: %d\n', ...
    nnz(hard),nnz(soft),nnz(classCode==2));
fprintf('Direct soft/hard grain contacts: %d | summed interface length: %.3f um\n', ...
    nContact,contactLen);
fprintf('All SF are geometric. Basal normal stress uses uniform sigma %.1f MPa.\n',sigma);

% FIG. 22: categorical grain map; highlighted RED alpha-alpha soft/hard contact.
try
    f=figure('Name','22 Soft-Hard Alpha Grain Interaction','Color','white');
    plot(gs,classCode,'micronbar','on'); hold on;
    colormap(gca,[0.25 0.47 0.86;0.69 0.69 0.69;0.96 0.70 0.16]);
    setColorRange([0.5 3.5]);
    cb=mtexColorbar; try,cb.Ticks=[1 2 3];cb.TickLabels=classNames;catch,end
    set(gca,'Color',[0.12 0.12 0.12]);
    plot(gs.boundary,'lineColor',[0.23 0.23 0.23],'lineWidth',0.35);
    if any(segSH), plot(B(segSH),'lineColor',[0.9 0.05 0.1],'lineWidth',2.4);end
    hold off; axis tight;
    title(sprintf('22 Soft-hard alpha grain contacts | [%g %g %g] | %d interfaces',...
        xyz,nContact));
    softHardExport(f,folder,'22_SoftHard_GrainInteraction');
catch ME
    warning('Figure 22 failed: %s',ME.message);
end

% FIG. 23: contrast on each genuine alpha-alpha boundary segment.
try
    f=figure('Name','23 Easy-a Schmid Contrast Across Grain Boundaries','Color','white');
    plot(gs,zeros(n,1),'micronbar','on'); hold on;
    colormap(gca,gray(256)); setColorRange([0 1]);
    if ~isempty(segLen)
        plot(B,segContrast,'lineWidth',2.2);
        colormap(gca,parula(256));setColorRange([0 0.5]);
        mtexColorbar('title','Delta m easy-<a> (absolute)');
        if any(segSH)
            plot(B(segSH),'lineColor',[1 0 0],'lineWidth',3.1);
        end
    end
    hold off;axis tight;set(gca,'Color',[0.12 0.12 0.12]);
    title('23 Alpha-alpha boundary contrast in max(basal,prismatic) SF | red = soft-hard');
    softHardExport(f,folder,'23_SchmidContrast_AlphaGrainBoundaries');
catch ME
    warning('Figure 23 failed: %s',ME.message);
end

% FIG. 24: basal normal-stress fraction per alpha GRAIN (c-axis driven).
try
    f=figure('Name','24 Basal Normal Stress Ratio by Alpha Grain','Color','white');
    plot(gs,sigN,'micronbar','on');hold on;
    plot(gs.boundary,'lineColor',[0.25 0.25 0.25],'lineWidth',0.4);
    if any(segSH),plot(B(segSH),'lineColor',[1 0 0],'lineWidth',2.5);end
    hold off; colormap(gca,parula(256));setColorRange([0 1]);
    mtexColorbar('title','sigma_n(basal)/sigma = cos^2(theta)');
    title('24 Basal-plane normal stress fraction | uniform tensile stress assumption');
    set(gca,'Color',[0.12 0.12 0.12]);
    softHardExport(f,folder,'24_BasalNormalStress_AlphaGrains');
catch ME
    warning('Figure 24 failed: %s',ME.message);
end

OUT=struct('LoadAxisXYZ',xyz,'GrainAngleDeg',opt.grainAngleDeg,...
    'MinGrainPixels',opt.minGrainPixels,...
    'NGrains',n,'NHard',nnz(hard),'NSoft',nnz(soft),'NIntermediate',nnz(classCode==2),...
    'NSoftHardContactPairs',nContact,'SoftHardContactLength_um',contactLen,...
    'GrainData',Tgr,'NeighborPairData',Tpair);
save(fullfile(folder,'22_SoftHard_Results.mat'),'OUT','opt','-v7.3');
end

function softHardExport(fig,folder,base)
try
    exportgraphics(fig,fullfile(folder,[base '.png']),'Resolution',260);
catch ME
    warning('PNG export problem: %s',ME.message);
end
try
    savefig(fig,fullfile(folder,[base '.fig']));
catch ME
    warning('FIG export problem: %s',ME.message);
end
end
