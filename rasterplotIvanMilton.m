function rasterplotIvanMilton(TrialData,varargin)
%
%  rasterplotgid(TrialData,OPTIONS)
%
%  OPTIONS:  
%     'SelectNeuron'  'spike11'
%     'Task_IDnumber' [1 2]
%     'SortTrials'    'yes'
%     'quickplot'     'yes'
%     'xlim'          [-1 2]
%
%  TrialData can be either a 'filename' or a trial structure.

% Set default parameter values or get them from VARARGIN.
% Will sort based on coherence, motion dir, and response (corr or wrong)
SortTrials = getArgumentValue('sorttrials','yes',varargin{:});

% Plots dots instead of lines.
quickplot  = getArgumentValue('quickplot','yes',varargin{:});
xlims       = getArgumentValue('xlim',[-inf inf],varargin{:});% Set xlim.
markersOn  = getArgumentValue('markersOn',true,varargin{:});  % Event markers

if length(TrialData)<1
    disp('No TrialData')
    return
end

togglefig('rasterplot'),clf, drawnow,pause(.1)
% subplot(3,2,[1,3,5])
% https://la.mathworks.com/matlabcentral/answers/348450-return-figure-handle-from-a-function
% Load the file if TrialData is a filename string.
if ischar(TrialData)
   t=load(TrialData);  % Load the trial data.
   trial = t.trials;

   Fields = fields(trial);                % Get the field names.
   % Display the SpikeCodes of the neurons recorded in this file.
   disp('The file contains the following spikeCodes: ') 
   indices = strmatch('spike',Fields);
   disp(Fields(indices))
   % Thsese inputs arguments are mandatory when reading data from a file. 
   % Get this for the figure title string
   neuronSelected = getArgumentValue('selectneuron',[],varargin{:});
   % Get this for the figure title string
   isRhythmTask  = getArgumentValue('rhythmTask',[],varargin{:});
   % Verify that the 'align' input is present.
   alignEvent = getArgumentValue('align',[],varargin);
%   togglefig([TrialData '    task ' num2str(isRhythmTask) ', ' neuronSelected])
%   set(gca,'position',[.13 .03 .86 .96])
else
   trial = TrialData;
end

% Select the trials based on the argument,[parameter] pairs in the input.
%trial = selectTrials(trial,varargin{:});

% Las columnas de la variable conducta son:
% (1) Tarea, 0:timing, 1:centerOut
% (2) esControl, 1:esControl, 0:noEsControl
% (3) ladoInicio, 1:izquierda, 2:derecha
% (4) duracion: 2:500ms, 3:750ms, 4:1000ms
% (5) numeroSamples
% (6) numero intervalos de memoria
% (7) esHit
% asi los ordena ordenaotto:[7 3 6 4 2 1]

sortmatrix = [ [trial.tGoCue]' [trial.cEsHit]' [trial.cIniciaIzq]' ...
               [trial.cNumeroMemoria]' [trial.cDurIntervalo]' ...
               [trial.cEsControl]' [trial.cEsTiming]' ];

sortmatrix = fliplr(sortmatrix);
          
sorted = 1:length(trial);
if strcmpi(SortTrials,'yes') % Sort the trials
   [~, sorted] = sortrows(sortmatrix);
end

%sortmatrix = flipud(sortmatrix);
sorted = flipud(sorted);

% Initialize the spikes and markers matrices.
spiketimes   = cell(length(trial),1);
markertimes1 = cell(length(trial),1);
markertimes2 = cell(length(trial),1);

% Loop through each trial to get the spikes and markers.
for k = 1:length(trial)
  
       
   ktrial = trial(sorted(k));
%    alignTime = ktrial.tGoCue;
   alignTime = 0;
   

   spiketimes{k} = (ktrial.spikeTimes(ktrial.spikeTimes>=xlims(1) & ...
                    ktrial.spikeTimes<=xlims(2)) ) - alignTime;

   if ktrial.cIniciaIzq==1 % Stim init left
%       markertimes1U   = [ktrial.tMuestraIzq(:); ...
%                          ktrial.tMuestraDer(:)]' - alignTime;
       markertimes1U   = ktrial.tMuestra(:) - alignTime;
       markertimes1{k} = markertimes1U(markertimes1U>=xlims(1) & ...
                                       markertimes1U<=xlims(2));
%       markertimes2U   = [ktrial.tMemoriaDer(:); ktrial.tMemoriaIzq(:); ...
%                          ktrial.tGoCue(:) ]' - alignTime;
       markertimes2U   = [ktrial.tMemoria(:); ktrial.tGoCue(:) ]' - alignTime;
       markertimes2{k} = markertimes2U(markertimes2U>=xlims(1) &  ...
                                       markertimes2U<=xlims(2));
   else
%       markertimes3U   = [ktrial.tMuestraIzq(:); ...
%                          ktrial.tMuestraDer(:)]' - alignTime;
       markertimes3U   = ktrial.tMuestra(:) - alignTime;
       markertimes3{k} = markertimes3U(markertimes3U>=xlims(1) & ...
                                       markertimes3U<=xlims(2));
%       markertimes4U   = [ktrial.tMemoriaDer(:); ktrial.tMemoriaIzq(:); ...
%                          ktrial.tGoCue(:) ]' - alignTime;
       markertimes4U   = [ktrial.tMemoria(:); ktrial.tGoCue(:) ]' - alignTime;
       markertimes4{k} = markertimes4U(markertimes4U>=xlims(1) & ...
                                       markertimes4U<=xlims(2));
   end
   
   if ktrial.cEsHit==0 % errors
       markertimes7{k} = 0;
   else
       markertimes7{k} = [];
   end

   markertimes5U   = ktrial.tIniMov - alignTime;
   markertimes5{k} = markertimes5U(markertimes5U>=xlims(1) &  ...
                     markertimes5U<=xlims(2));

   markertimes6U   = ktrial.tFinMov - alignTime;
   markertimes6{k} = markertimes6U(markertimes6U>=xlims(1) &  ...
                     markertimes6U<=xlims(2));
end

if strcmp(quickplot,'yes')
   linestyle = 'none'; marker    = '.';
else
   linestyle = '-'; marker    = 'none';
end

[spikesx, spikesy]   = rasterplot(spiketimes  ,'QuickPlot',quickplot);

% Plot the markers
if markersOn
    [markers1x, markers1y] = rasterplot(markertimes1,'QuickPlot',quickplot);
    [markers2x, markers2y] = rasterplot(markertimes2,'QuickPlot',quickplot);
    [markers3x, markers3y] = rasterplot(markertimes3,'QuickPlot',quickplot);
    [markers4x, markers4y] = rasterplot(markertimes4,'QuickPlot',quickplot);
    [markers5x, markers5y] = rasterplot(markertimes5,'QuickPlot',quickplot);
    [markers6x, markers6y] = rasterplot(markertimes6,'QuickPlot',quickplot);
    [markers7x, markers7y] = rasterplot(markertimes7,'QuickPlot',quickplot);
    line(markers2x,markers2y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[56,77,82]/255), % memoria izquierda
    line(markers1x,markers1y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[131,6,137]/255)  %entrainment izquierda
    line(markers4x,markers4y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[95,163,103]/255) % memoria derecha
    line(markers3x,markers3y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[162,136,227]/255) % entrainment derecha
    line(markers5x,markers5y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[64,207,237]/255)    % gocue
    line(markers6x,markers6y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color',[3,157,255]/255)    %reward
    line(markers7x,markers7y,'LineStyle',linestyle,'marker',marker, ...
        'markersize',15,'linewidth',3,'color','r')
end

% Set figure properties.
% set(gcf,'units','normalized')
% set(gcf,'position',[.5 0 .5 .99])
% set(gca,'ylim',[0 length(trial)+1],'ytick',1:1:length(trial),'xlim',xlim,...
%         'yticklabel',num2str(round(sortmatrix(1:1:end,:)*10)/10),'fontsize',6, 'tickdir','out'); box off
% set(gca,'ylim',[0 length(trial)+1],'ytick',[],'xlim',xlim,...
%         'fontsize',6, 'tickdir','out','yticklabel',[]); box off
%set(gca,'ylim',[0 length(trial)+1],'fontsize',6, 'tickdir','out','ytick',[[1 5:5:length(trial)+1]]); box off
% Sin ejes
% set(gca,'ylim',[0 length(trial)+1],'ytick',[],...
%         'fontsize',6, 'tickdir','out','yticklabel',[]); box off

% Plot the spikes
line(spikesx, spikesy,   'LineStyle',linestyle,'marker',marker,'markersize', 10,'linewidth',1,'color','k')

% if ktrial.spikeTimes


set(gcf,'color','w')
set(gca,'color','w')
set(gcf,'position',[   494   202   777   503])
set(gca,'FontSize',18,'fontname','times')
ylim([0 k]), xlim(xlims);
xlabel('time (sec)')
ylabel('trials')

box off, grid off
shg

x = [0,1];
y = [-1,-1];
colores = {[131,6,137]/255  ,  [162,136,227]/255  ,  [56,77,82]/255  ,  [95,163,103]/255,...
    [64,207,237]/255,[3,157,255]/255};
h = gobjects(length(colores),1);
for f = 1:length(colores)
    h(f) = line(x,y,'color',colores{f},'linewidth',3);
    hold on
end
legend(h,{'entrainment L','entrainment R','memory L','memory R','InitMov','FinMov'},'location','northeast')

h = togglefig('rasterplot');
% raster  = ancestor(h, 'figure');
set(gcf, 'Units','normalized','Position',[0.5 0.5 0.5 0.5]);

drawnow