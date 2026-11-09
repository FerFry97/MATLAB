%% Ti64 CHANNEL 5 vs MTEX: phase-complete measured texture validation
% Stand-alone. Does NOT reconstruct parent beta and does NOT discard any
% indexed alpha/beta pixel. Tested by code review only (MATLAB unavailable).
% Two approaches: (A) MTEX ODF -> PF and (B) direct specimen pole density.
% B is not claimed to duplicate the proprietary CHANNEL 5 density estimator.
% Use the exported CTF for the SAME Channel 5 selection/subset.
clearvars; close all; clc;

%% SETTINGS -- EDIT THESE IF NEEDED
mtexRoot = 'C:\mtex-6.1.1'; % startup_mtex.m directory
halfwidthDeg = 10;
runDirectPF = true;       % direct measured-poles density; can be memory-intensive
runSensitivity = false;   % additional halfwidths [5,7.5,10,12.5,15] deg
% Coordinate treatment of CTF: compare both only if needed and keep them
% in separate runs. Never rotate Euler data just to match intensity numbers.
ctfConversion = 'convertEuler2SpatialReferenceFrame';
channel5AlphaBasalMUD = 19.35; % Sample 8 0.15 um Subset2, clustering ON at 5 deg
channel5AlphaBasalNoClusterMUD = 18.88;

if exist('EBSD','class') ~= 8
    startupFile = fullfile(mtexRoot,'startup_mtex.m');
    if exist(startupFile,'file') ~= 2
        error('MTEX not on path. Set mtexRoot to your MTEX installation root.');
    end
    run(startupFile);
end
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

[fn,fp] = uigetfile('*.ctf','Choose the exact CTF exported from CHANNEL 5');
if isequal(fn,0), return; end
f = fullfile(fp,fn);
[~,shortName] = fileparts(f);
out = fullfile(fp,[shortName '_Channel5_MTEX_Texture_Validation']);
if ~exist(out,'dir'), mkdir(out); end

switch ctfConversion
    case 'convertEuler2SpatialReferenceFrame'
        E = EBSD.load(f,'interface','ctf','convertEuler2SpatialReferenceFrame');
    case 'convertSpatial2EulerReferenceFrame'
        E = EBSD.load(f,'interface','ctf','convertSpatial2EulerReferenceFrame');
    otherwise
        error('Unsupported ctfConversion setting.');
end
alpha = E('Ti-Hex');
beta  = E('Titanium cubic');
fprintf('CTF: %s | total %d | alpha %d | measured beta %d | unindexed %d\n', ...
  fn,length(E),length(alpha),length(beta),length(E)-length(E('indexed')));
fprintf('Halfwidth = %.2f deg; reference mode = %s\n',halfwidthDeg,ctfConversion);
if isempty(alpha), error('No indexed Ti-Hex data in selected CTF.'); end

% phase-specific plane families (hex 4-index and cubic 3-index).
hA = [Miller(0,0,0,2,alpha.CS), Miller(1,0,-1,0,alpha.CS), ...
      Miller(1,1,-2,0,alpha.CS)];
namesA = {'0002','10-10','11-20'};
if ~isempty(beta)
    hB = [Miller(1,1,0,beta.CS), Miller(2,0,0,beta.CS), ...
          Miller(2,1,1,beta.CS)];
    namesB = {'110','200','211'};
end

% Long format avoids concatenation dimension errors; maxima, not colorbar top.
T = table('Size',[0,7], ...
 'VariableTypes',{'string','string','string','double','double','double','double'}, ...
 'VariableNames',{'Phase','Method','PF','PixelCount','Halfwidth_deg','Max_MUD_MRD','RatioToC5Basal'});
T = [T; runPhase(alpha,hA,namesA,'alpha',halfwidthDeg,runDirectPF,out,channel5AlphaBasalMUD)];
if ~isempty(beta)
    T = [T; runPhase(beta,hB,namesB,'measured_beta',halfwidthDeg,runDirectPF,out,NaN)];
end
writetable(T,fullfile(out,'AllMeasuredPF_Peaks.csv'));
disp(T);

if runSensitivity
    Tsen = table();
    for hw = [5 7.5 10 12.5 15]
        % Only ODF estimation sensitivity; not a CHANNEL 5 cluster-size scan.
        tmp = runPhase(alpha,hA(1),namesA(1),'alpha_sensitivity',hw,false,out,channel5AlphaBasalMUD);
        Tsen = [Tsen;tmp]; %#ok<AGROW>
    end
    writetable(Tsen,fullfile(out,'AlphaBasal_HalfwidthSensitivity.csv'));
end

fprintf('\nCHANNEL 5 reference: alpha basal 19.35 MUD (cluster 5 deg); 18.88 MUD (no cluster).\n');
fprintf('MTEX direct density != proprietary CHANNEL 5 algorithm. Compare shapes, peak positions, normalization.\n');
fprintf('Output: %s\n',out);

%% FUNCTIONS
function T = runPhase(data,h,names,phase,hw,direct,folder,c5max)
  N = length(data);
  ori = data.orientations;
  odf = calcDensity(ori,'halfwidth',hw*degree,'silent');
  T = table('Size',[0,7], ...
   'VariableTypes',{'string','string','string','double','double','double','double'}, ...
   'VariableNames',{'Phase','Method','PF','PixelCount','Halfwidth_deg','Max_MUD_MRD','RatioToC5Basal'});
  for k=1:length(h)
    label = string(names{k});
    % ODF projected to pole figure, using all pixels in this phase
    pf = calcPDF(odf,h(k));
    v = max(pf);
    T = [T; makeRow(phase,'ODF_to_PF',label,N,hw,double(v),c5max)]; %#ok<AGROW>
    fh = figure('Visible','off','Color','white','Name',sprintf('%s %s ODF',phase,names{k}));
    try
        plotPDF(odf,h(k),'antipodal','contourf','upper','projection','equalarea');
    catch
        plotPDF(odf,h(k),'antipodal','contourf');
    end
    mtexColorbar('title','MRD');
    title(sprintf('%s {%s} | ODF->PF | N=%d | hw=%.1f deg | max=%.3f',phase,names{k},N,hw,double(v)), ...
     'Interpreter','none');
    exportgraphics(fh,fullfile(folder,sprintf('%s_%s_ODF_PF.png',phase,names{k})),'Resolution',220);
    close(fh);

    if direct
      try
        % Symmetrically equivalent poles for each measured orientation.
        % Critical: mean-one normalized density of equivalent directions.
        hk = symmetrise(h(k),'unique','antipodal');
        poles = ori * hk;
        poles.antipodal = true;
        pdfDirect = calcDensity(poles,'halfwidth',hw*degree,'antipodal');
        vmax = max(pdfDirect);
        T = [T;makeRow(phase,'Direct_pole_KDE',label,N,hw,double(vmax),c5max)]; %#ok<AGROW>
        fh=figure('Visible','off','Color','white','Name',sprintf('%s %s Direct',phase,names{k}));
        plot(pdfDirect,'antipodal','contourf','upper');
        mtexColorbar('title','MUD');
        title(sprintf('%s {%s} | Direct pole density | N=%d | hw=%.1f deg | max=%.3f', ...
          phase,names{k},N,hw,double(vmax)),'Interpreter','none');
        exportgraphics(fh,fullfile(folder,sprintf('%s_%s_Direct_PF.png',phase,names{k})),'Resolution',220);
        close(fh);
      catch ME
        warning('Direct density failed for %s %s: %s',phase,names{k},ME.message);
      end
    end
  end
end
function row = makeRow(phase,method,pf,N,hw,v,c5)
  ratio=NaN;
  if isfinite(c5) && pf == "0002", ratio=v/c5; end
  row = {string(phase),string(method),string(pf),double(N),double(hw),double(v),double(ratio)};
  row = cell2table(row,'VariableNames', ...
    {'Phase','Method','PF','PixelCount','Halfwidth_deg','Max_MUD_MRD','RatioToC5Basal'});
end
