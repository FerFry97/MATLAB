%% Ti-6Al-4V: pole-figure intensity validation, MTEX 6.1.1
% Stages covered: 1 (data audit), 2 (reproduce MTEX maxima), 3 (halfwidth
% sensitivity), 4 (normalisation checks), 5 (peak location).
% NOT TESTED in MATLAB (written without MATLAB access). The same algorithm was
% run in Python on SAMPLE015.ctf and reproduced your MTEX maxima within ~0.2%.
% Nothing here modifies orientations; all numbers come from the raw CTF.
%
% Output CSVs go to outputFolder. Channel 5 values can be typed in section 1.

clear; close all; clc;

%% 0. MTEX
if exist('EBSD','class') ~= 8
    cand = {'C:\mtex-6.1.1\startup_mtex.m','C:\mtex-6.1.1\startup.m'};
    ok = false;
    for k = 1:numel(cand)
        if exist(cand{k},'file') == 2, run(cand{k}); ok = true; break, end
    end
    if ~ok, error('MTEX not found. Edit the candidate paths.'); end
end
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
setMTEXpref('defaultColorMap',WhiteJetColorMap);

%% 1. Options
alphaPhase = 'Ti-Hex';  betaPhase = 'Titanium cubic';
expectedAlphaN = 134025;  expectedBetaN = 860;     % from the CTF audit
halfwidths_deg = [2.5 5 7.5 10 12.5 15 20];        % Stage 3
refHalfwidth_deg = 10;                             % Stage 2
% Reference maxima to reproduce (your AllMeasuredPF_Peaks.csv values)
ref.alpha = [7.648542 3.272208 3.465820];          % {0002} {10-10} {11-20}
ref.beta  = [2.300323 4.141814 1.997100];          % {110} {200} {211}
ch5_alpha0002 = [18.88 19.35];                     % Channel 5, no clustering / clustering
gridRes = 1*degree;                                % evaluation grid for ODF-based PF

[fileName,filePath] = uigetfile({'*.ctf','Oxford CTF (*.ctf)'},'Select the CTF');
if isequal(fileName,0), error('No file selected.'); end
fname = fullfile(filePath,fileName);
outputFolder = fullfile(filePath,'PF_validation_MTEX611');
if ~isfolder(outputFolder), mkdir(outputFolder); end

%% 2. Stage 1: data audit
ebsd = EBSD.load(fname,'interface','ctf','convertEuler2SpatialReferenceFrame');
alphaE = ebsd(alphaPhase);  betaE = ebsd(betaPhase);
CSa = alphaE.CS;  CSb = betaE.CS;
fprintf('alpha N = %d (expected %d) | beta N = %d (expected %d)\n', ...
    length(alphaE),expectedAlphaN,length(betaE),expectedBetaN);
if length(alphaE) ~= expectedAlphaN || length(betaE) ~= expectedBetaN
    warning('Pixel counts differ from the audited CTF: you are not using the same data set.');
end
disp(CSa); disp(CSb);   % check 6/mmm and m-3m, and the hexagonal axis convention
try, disp(CSa.alignment); catch, end   % hexagonal axis convention (a vs a*)

%% 3. Pole-figure families
hA = {Miller(0,0,0,2,CSa),Miller(1,0,-1,0,CSa),Miller(1,1,-2,0,CSa)};
nA = {'alpha {0002}','alpha {10-10}','alpha {11-20}'};
hB = {Miller(1,1,0,CSb),Miller(2,0,0,CSb),Miller(2,1,1,CSb)};
nB = {'beta {110}','beta {200}','beta {211}'};
% Each field is wrapped in an extra {} so that "sets" is ONE struct whose
% fields are cell arrays (otherwise struct() builds a 1x6 struct array).
oriList = [repmat({alphaE.orientations},1,3) repmat({betaE.orientations},1,3)];
sets = struct('name',{[nA nB]},'h',{[hA hB]},'ori',{oriList});

r = equispacedS2Grid('upper','resolution',gridRes);   % evaluation directions

%% 4. Stage 2 + 3: ODF-based maxima (MTEX) over halfwidths
rows = {};
for iw = 1:numel(halfwidths_deg)
    hw = halfwidths_deg(iw)*degree;
    for s = 1:numel(sets.name)
        odf = calcDensity(sets.ori{s},'halfwidth',hw,'silent');
        pdf = odf.calcPDF(sets.h{s},r);
        [m,id] = max(pdf(:));
        rp = r(id);
        rows(end+1,:) = {sets.name{s},halfwidths_deg(iw),'MTEX ODF->PF',m, ...
            rp.theta/degree,mod(rp.rho/degree,360)}; %#ok<AGROW>
    end
    fprintf('halfwidth %.1f deg done\n',halfwidths_deg(iw));
end

%% 5. Direct (ODF-free) estimators: cone and von Mises-Fisher on the sphere
% Independent check of normalisation: density is relative to a uniform
% distribution on the sphere (uniform = 1), antipodal (axial) poles.
for iw = 1:numel(halfwidths_deg)
    hw = halfwidths_deg(iw)*degree;
    for s = 1:numel(sets.name)
        X = poleXYZ(sets.ori{s},sets.h{s});
        [mc,c1] = maxDirect(X,'cone',hw);
        [mv,c2] = maxDirect(X,'vmf',hw);
        rows(end+1,:) = {sets.name{s},halfwidths_deg(iw),'direct cone',mc, ...
            acosd(abs(c1(3))),mod(atan2d(c1(2)*sign(c1(3)),c1(1)*sign(c1(3))),360)}; %#ok<AGROW>
        rows(end+1,:) = {sets.name{s},halfwidths_deg(iw),'direct vMF (HWHM)',mv, ...
            acosd(abs(c2(3))),mod(atan2d(c2(2)*sign(c2(3)),c2(1)*sign(c2(3))),360)}; %#ok<AGROW>
    end
end

T = cell2table(rows,'VariableNames',{'PoleFigure','Halfwidth_deg','Method','MaxIntensity', ...
    'PeakPolar_deg','PeakAzimuth_deg'});
writetable(T,fullfile(outputFolder,'PF_maxima_all_methods.csv'));

%% 6. Stage 2 report: reproduction of the reference
fprintf('\n--- Reproduction at halfwidth %g deg ---\n',refHalfwidth_deg);
refv = [ref.alpha ref.beta];
for s = 1:numel(sets.name)
    k = strcmp(T.PoleFigure,sets.name{s}) & T.Halfwidth_deg == refHalfwidth_deg & ...
        strcmp(T.Method,'MTEX ODF->PF');
    fprintf('%-14s ref %.4f | now %.4f | diff %.2f%%\n',sets.name{s},refv(s), ...
        T.MaxIntensity(k),100*(T.MaxIntensity(k)/refv(s)-1));
end

%% 7. Stage 3 figure: maximum intensity vs halfwidth
try
    figure('Color','white','Name','Max PF intensity vs halfwidth');
    meth = {'MTEX ODF->PF','direct cone','direct vMF (HWHM)'};
    for p = 1:2
        subplot(1,2,p); hold on;
        for m = 1:numel(meth)
            k = strcmp(T.PoleFigure,sets.name{(p-1)*3+1+(p==2)}) & strcmp(T.Method,meth{m});
            loglog(T.Halfwidth_deg(k),T.MaxIntensity(k),'o-','DisplayName',meth{m});
        end
        set(gca,'XScale','log','YScale','log'); grid on;
        xlabel('halfwidth [deg]'); ylabel('max intensity [mrd]');
        title(sets.name{(p-1)*3+1+(p==2)}); legend('Location','northeast');
    end
    if exist('exportgraphics','file') == 2
        exportgraphics(gcf,fullfile(outputFolder,'Max_vs_halfwidth.png'),'Resolution',300);
    end
catch ME
    warning('Figure failed: %s',ME.message);
end

%% 8. Stage 4: normalisation checks
odf = calcDensity(alphaE.orientations,'halfwidth',refHalfwidth_deg*degree,'silent');
rr = equispacedS2Grid('resolution',2*degree);          % whole sphere
pdfAll = odf.calcPDF(hA{1},rr);
fprintf('{0002} mean over the full sphere (should be 1): %.4f\n',mean(pdfAll(:)));

%% 9. Channel 5 comparison for {0002}
k = strcmp(T.PoleFigure,'alpha {0002}') & strcmp(T.Method,'MTEX ODF->PF');
hw = T.Halfwidth_deg(k);  v = T.MaxIntensity(k);
for t = ch5_alpha0002
    heq = exp(interp1(log(v),log(hw),log(t)));
    fprintf('Channel 5 max %.2f corresponds to an MTEX halfwidth of ~%.2f deg (interpolated)\n',t,heq);
end
fprintf('This is a DESCRIPTIVE equivalence, not proof of how Channel 5 defines its halfwidth.\n');
fprintf('\nDONE. Output folder: %s\n',outputFolder);

%% ===================== local functions =====================
function X = poleXYZ(ori,h)
% sample-frame unit vectors of all symmetrically equivalent poles, axial
hs = symmetrise(h);  hs = hs(:);
X = zeros(0,3);
for j = 1:numel(hs)
    v = ori * hs(j);
    X = [X; v.xyz]; %#ok<AGROW>
end
X = X ./ vecnorm(X,2,2);
end

function [mx,best] = maxDirect(X,kind,hw)
% Maximum of an ODF-free axial pole-density estimate (uniform = 1)
n = size(X,1);
if n > 40000   % thin the data for speed, keeping the density scale
    idx = randperm(n,40000);  Xs = X(idx,:);
else
    Xs = X;
end
switch kind
    case 'cone'
        c0 = cos(hw);
        dens = @(C) sum(abs(C*Xs') >= c0,2) / size(Xs,1) / (1-c0);
    case 'vmf'
        k = log(2)/(1-cos(hw));
        a = k/(1-exp(-2*k));    % axial (antipodal) normalisation: uniform = 1
        dens = @(C) a * mean(exp(k*(abs(C*Xs')-1)) + exp(k*(-abs(C*Xs')-1)),2);
end
cand = Xs(randperm(size(Xs,1),min(1500,size(Xs,1))),:);
d = zeros(size(cand,1),1);
for i = 1:200:size(cand,1)
    j = i:min(i+199,size(cand,1));  d(j) = dens(cand(j,:));
end
[~,o] = sort(d,'descend');  mx = 0;  best = cand(o(1),:);
for q = 1:3
    c = cand(o(q),:);  v = d(o(q));
    for sdeg = [2 1 0.5 0.25 0.1 0.05]
        R = c + deg2rad(sdeg)*randn(150,3);
        R = [c; R ./ vecnorm(R,2,2)];
        dd = dens(R);  [v,jj] = max(dd);  c = R(jj,:);
    end
    if v > mx, mx = v; best = c; end
end
end
