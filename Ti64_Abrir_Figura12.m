function Ti64_Abrir_Figura12()
% Displays the saved six-PF alpha / reconstructed-beta comparison in MATLAB.
% Run from Command Window: Ti64_Abrir_Figura12
% No MTEX or EBSD reprocessing is needed for viewing an existing PNG/FIG.
%
% MATLAB Editor opens .m source code. Actual plots appear in a Figure
% window or in the desktop's docked Figures panel.

    [fileName,folder] = uigetfile( ...
        {'*.fig;*.png','Figure 12 files (*.fig, *.png)'; ...
         '*.fig','Editable MATLAB figure (*.fig)'; ...
         '*.png','Saved PNG (*.png)'}, ...
        'Select 12_PF_Alpha_vs_Beta_MatchedPixels');

    if isequal(fileName,0)
        fprintf('Canceled: no file selected.\n');
        return;
    end

    selected = fullfile(folder,fileName);
    [~,stem,ext] = fileparts(selected);
    pngFile = fullfile(folder,[stem '.png']);
    figFile = fullfile(folder,[stem '.fig']);
    if ~contains(fileName,'12_PF_Alpha_vs_Beta_MatchedPixels')
        warning('You selected %s, not the expected Figure 12.',fileName);
    end

    h = [];
    % Try the editable FIG before displaying a flat image.
    if exist(figFile,'file') == 2
        try
            h = openfig(figFile,'new','visible');
            % A correct editable 6-PF plot should contain >= 6 axes.
            numAxes = numel(findall(h,'Type','axes'));
            if numAxes < 6
                warning(['Saved FIG contains only %d axes, expected 6. ', ...
                    'Falling back to the PNG if possible.'],numAxes);
                close(h);
                h = [];
            else
                fprintf('Opened editable FIG with %d axes.\n',numAxes);
            end
        catch ME
            warning('Could not open editable .fig: %s',ME.message);
            h = [];
        end
    end

    if isempty(h)
        if exist(pngFile,'file') ~= 2
            error(['Neither an editable six-PF FIG nor a PNG could be opened.\n' ...
                   'Expected files in: %s'],folder);
        end
        rgb = imread(pngFile);
        h = figure('Name','12 | Alpha vs reconstructed beta | six PF (PNG)', ...
                   'NumberTitle','off','Color','white','Visible','on');
        ax = axes('Parent',h);
        image(ax,rgb);
        axis(ax,'image'); axis(ax,'off');
        title(ax,'Saved six pole figures (image, not editable PF data)', ...
              'Interpreter','none');
        fprintf('Opened PNG. Individual PFs are not editable in an image.\n');
    end

    set(h,'Visible','on');
    try
        set(h,'WindowState','maximized');
    catch
        % Older MATLAB versions may not have WindowState.
    end
    figure(h);
    drawnow;
    shg;

    fprintf('\nIf you want the figure DOCKED in MATLAB desktop, run:\n');
    fprintf('    set(gcf,''WindowStyle'',''docked'');\n\n');
    fprintf('Figure loaded from: %s\n',folder);
end
