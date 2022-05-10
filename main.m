
%                                   Main 
%                  »»»»»»»»»»»  ver May 2022  «««««««««««

% ░░░░░░░░░░░░░░░░░░░░░░░░░░░   ~ code ~   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
% Contiene : 
%   - main paths en diferentes pcs
%   - loads paleta de colores en 'coolors.m'
%   - loads plot presets (font size, etc) 
%   - loads indexes for data analysis (of the RhythmTask structure)
% ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░

% ╔════════════════════════════════════════════════════════════════════════════╗
% ║ TODO:                                                                      ║ 
%   selecionar lista con codigo de refinaListaAreas
%   lista   = refinaListaAreas(recArea,'listaForPSD',false,'tipoN',1);
% ╚════════════════════════════════════════════════════════════════════════════╝


% % Path
if exist('/Users/laboratorioB05/') % Mac Lab
    
    disp('accediendo a directorio LafuenteLab')
    
    mainpath = '/Users/laboratorioB05/Documents/phd/MATLAB/scriptsMILI/';
    figpath = '/Users/laboratorioB05/Documents/phd/MATLAB/scriptsMILI/figuras/';
    datapath =...
        '/Users/laboratorioB05/Documents/phd/MONOS/DATOS_IVAN_18_10_2019/EstructurasIvan/matFilesHippocampus/';
    addpath(genpath('/Users/laboratorioB05/Documents/phd/MATLAB/scriptsMILI/funMili/'))
    
    addpath(genpath(mainpath));
    addpath(genpath(datapath));
    addpath(genpath('/Users/laboratorioB05/Documents/phd/MATLAB/descargas_mathworks')) 
    
    if exist('/Volumes/VIRUS')
        disp('loading external hard drive')
        addpath(genpath('/Volumes/VIRUS/datosIvanMilton/estructurasOtto'))
    end    
    
elseif exist('/Users/administrador/Music/DOCTORADO/') % MacBookPro
    
    disp('accediendo a directorio laptop')
    
    mainpath = '/Users/administrador/Music/DOCTORADO/MATLAB/scriptsMILI_MacBookPro/';
    figpath = '/Users/administrador/Music/DOCTORADO/MATLAB/scriptsMILI_MacBookPro/figuras/';
    datapath =...
        '/Users/administrador/Music/DOCTORADO/MONOS/IVAN/DATOS_IVAN_18_10_2019/EstructurasIvan/matFilesHippocampus/';
   
    addpath(genpath(mainpath));
    addpath(genpath(datapath));
    addpath('/Users/administrador/Music/DOCTORADO/MATLAB/descargas_mathworks/')
    
    cd '/Users/administrador/Music/DOCTORADO/MATLAB/scriptsMILI_MacBookPro/'
    
    if exist('/Volumes/VIRUS')
     disp('loading external hard drive')
     addpath(genpath('/Volumes/VIRUS/datosIvanMilton/estructurasOtto/'))
    end
    
end

% load colors presets for plots -----------------------------------------
run('coolors.m')

% % Area selection ---------------------------------------------------------
brainAreas      = {'hipp', 'sma', 'v4', 'lip', 'mip', 'pfc', 'm1'};


if area == 1
    % HIPP   5, 65, 69, 85  
    neurona = 1;
    load('listaHippoParaMildred.mat')
%     load('listaHippoParaMildredHippo.mat')
%     load('listaHipocampo.mat')

elseif area == 2
    % SMA    1, 14, 15, 20, 22, 23
    neurona = 1;
    addpath('/Volumes/virus/datosIvanMilton/estructurasOtto')
    datapath =  '/Volumes/virus/datosIvanMilton/estructurasOtto/';
    load('listaParaMildredSMA.mat')  % lista de Victor

elseif area == 3
    % V4     9, 14, 16, 19, 28
    neurona = 1;
    addpath('/Volumes/virus/datosIvanMilton/estructurasOtto')
    datapath =  '/Volumes/virus/datosIvanMilton/estructurasOtto/';
    load('listaV4.mat')
%     load('neuronasV4.mat')

end


% load plot presets  ------------------------------------------------------------
IniciaLado     = [ 1 1 1 0 0 0 ];  % 1 inicia izquierda
durInt     = [ .5 .75  1 .5 .75  1 ];
Interv     = { 'izq500' , 'izq750' , 'izq1000' , 'der500' , 'der750' , 'der1000' };
names      = {'500ms','750ms','1000ms'};
subp1      = [ 1 3 5 2 4 6 ];  % types of trials
subp2      = [ 1 5 9 2 6 10 ]; % up plots 
subp3      = [ 3 7 11 4 8 12 ]; % down plots 
subp4      = [311, 312, 313, 311, 312, 313];
subp5      = [311, 311, 312, 312, 313, 313];
stimsLines = [ {500:500:500*8.5}, {750:750:750*8.5}, {1000:1000:1000*8.5},...
             {500:500:500*8.5}, {750:750:750*8.5}, {1000:1000:1000*8.5} ];
factor     = [ 0.5/0.5, 0.5/0.75, 0.5/1, 0.5/0.5, 0.5/0.75, 0.5/1, ];


% parameters for 'population_main.m' 
% p.Inicia =  [1   0   1  0  1 0]; % 1 inicia izquierda
% p.durInt =  [.5 .5 .75 .75 1 1];
% p.tFin   = (5.5+3) .* p.durInt; 
% p.names    = {'500ms','750ms','1000ms'};
% p.subplot  = [311 311 312 312 313 313];
% p.row      = [1 1 2 2 3 3]; 
% p.lineaLR    =  {'b-','r-','b-','r-','b-','r-'};


% Firing rate ------------------------------------------------------------
tStep        = 0.010; % Displacement of the FR window
minFR        = 0;

if area == 1 || area == 3
    timeConstant = 0.100;  % Time constant for the FR window 
    tFinal       = (5.5+3) .* durInt; 
    tStop      = durInt .* tFinal(6) * 1000; % for LFPS
    tAxis      = { 0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3),...
                0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3) };
    if exist('params_pca_wraped','var')
        tFinal     = repmat((5.5+3) .* .5,1,6); 
        tAxis    = { 0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3),...
                    0:tStep:tFinal(4), 0:tStep:tFinal(5), 0:tStep:tFinal(6) };
    end
%     timeConstant = 0.05;  % Time constant for the FR window 
elseif area == 2
%     timeConstant = 0.1;  
    tFinal       = repmat((3.5+3) .* .5,1,6);  
    tAxis      = { 0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3),...
                0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3) };
    if exist('params_pca_wraped','var')
        tFinal     = repmat((3.5+3) .* .5,1,6); 
        tAxis    = { 0:tStep:tFinal(1), 0:tStep:tFinal(2), 0:tStep:tFinal(3),...
                    0:tStep:tFinal(4), 0:tStep:tFinal(5), 0:tStep:tFinal(6) };
    end
end

% % LFP --------------------------------------------------------------------
% Fs          = 1000;
% bands       = { 'all', 'Delta', 'Theta', 'Beta', 'Gamma1', 'Gamma2' };
% evaluateFrec= [ {1:65} ,{1:4} , {4:12} , {12:32} , {32:55} , {65:100} ];





% load indexes ------------------------------------------------------------
% entrainment all elements
indice.ent = [{1:find(tAxis{1} == durInt(1)*3)},...
    {1:find(tAxis{2} == durInt(2)*3)},...
    {1:find(tAxis{3} == durInt(3)*3)},...
    {1:find(tAxis{1} == durInt(1)*3)},...
    {1:find(tAxis{2} == durInt(2)*3)},...
    {1:find(tAxis{3} == durInt(3)*3)}];

indice.ent_stims = [{0, find(tAxis{1} == durInt(1)*1), find(tAxis{1} == durInt(1)*2),...
    find(tAxis{1} == durInt(1)*3)},...
    {0, find(tAxis{2} == durInt(2)*1), find(tAxis{2} == durInt(2)*2),...
    find(tAxis{2} == durInt(2)*3)},...
    {0, find(tAxis{3} == durInt(3)*1), find(tAxis{3} == durInt(3)*2),...
    find(tAxis{3} == durInt(3)*3)},...
    ];
% entrainment by stims




% comandos generales  ----------------------------------------------------
% Comandos que se cargan al iniciar matlab
% todo: tener unidades normalizadas por default? 
set(0,'Defaultaxesxcolor',ones(1,3).*.1)
set(0,'Defaultaxesycolor',ones(1,3).*.1)
set(0,'Defaultaxeszcolor',ones(1,3).*.1)
set(0,'defaulttextcolor', ones(1,3).*.1)
set(0,'DefaultAxesXGrid','off','DefaultAxesYGrid','off')
set(0,'defaultfigureposition',[1 1 500 500]') 
set(0,'defaultAxesFontSize',27)
set(0,'defaultfigurecolor','w')
set(groot,'defaultaxesbox','off')
set(0,'DefaultaxesLineWidth', .9) 
set(groot, 'defaultAxesTickDir', 'out');
set(0,'DefaultaxesFontName', 'Arial') 
set(0,'DefaultlegendFontName', 'Arial')
format compact


% PCA analisis -----------------------------------------------------------
if exist('params_pca', 'var') 

        tFin_ent   = [3 .* [.5 .75 1], 3 .* [.5 .75 1]]; % tiempo final para Entrainment
        tAxis_ent  = { 0:tStep:tFin_ent(1), 0:tStep:tFin_ent(2), 0:tStep:tFin_ent(3),...
            0:tStep:tFin_ent(1), 0:tStep:tFin_ent(2), 0:tStep:tFin_ent(3)  };
        idx_tAxis_ent_stims = { find(ismember(tAxis_ent{1},0:0.5:tAxis_ent{1}(end))),...
            find(ismember(tAxis_ent{2},0:0.75:tAxis_ent{2}(end))),...
            find(ismember(tAxis_ent{3},0:1:tAxis_ent{3}(end))),...
            find(ismember(tAxis_ent{1},0:0.5:tAxis_ent{1}(end))),...
            find(ismember(tAxis_ent{2},0:0.75:tAxis_ent{2}(end))),...
            find(ismember(tAxis_ent{3},0:1:tAxis_ent{3}(end)))
            };
        tAxis_mem  = {tFin_ent(1):tStep:tFinal(1), tFin_ent(2):tStep:tFinal(2),...
            tFin_ent(3):tStep:tFinal(3), tFin_ent(1):tStep:tFinal(1),...
            tFin_ent(2):tStep:tFinal(2), tFin_ent(3):tStep:tFinal(3)};
        idx_tAxis_mem_stims = { find(ismember(tAxis_mem{1},0:0.5:tAxis_mem{1}(end))),...
            find(ismember(tAxis_mem{2},0:0.75:tAxis_mem{2}(end))),...
            find(ismember(tAxis_mem{3},0:1:tAxis_mem{3}(end))),...
            find(ismember(tAxis_mem{1},0:0.5:tAxis_mem{1}(end))),...
            find(ismember(tAxis_mem{2},0:0.75:tAxis_mem{2}(end))),...
            find(ismember(tAxis_mem{3},0:1:tAxis_mem{3}(end)))
            };
        s1 = length(tAxis{1}); s2 = length(tAxis{2}); s3 = length(tAxis{3});
        lims = [s1, s2, s3, s1, s2, s3];

        % total limites [inferior, superior] de cada frecuencia
        lim = [ 1                        , s1; 
                s1 + 1                   , s1 + s2 ;...
                s1 + s2 + 1              , s1 + s2 + s3 ;...
                s1 + s2 + s3 + 1         , (s1*2) + s2 + s3 ;...
                (s1*2) + s2 + s3 + 1     , (s1*2) + (s2*2) + s3 ;...
                (s1*2) + (s2*2) + s3 + 1 , (s1*2) + (s2*2) + (s3*2)
                ];

        % entrainment limites [inferior, superior]
        lim_ent = [ lim(1,1) + length(tAxis_ent{1}) - 1;...
                lim(2,1) + length(tAxis_ent{2}) - 1;...
                lim(3,1) + length(tAxis_ent{3}) - 1;...
                lim(4,1) + length(tAxis_ent{1}) - 1;...
                lim(5,1) + length(tAxis_ent{2}) - 1;...
                lim(6,1) + length(tAxis_ent{3}) - 1;...
                ];


        % % IDX CONTINUO
        idx_ent = [ lim(1,1) : lim_ent(1),...
                lim(2,1) : lim_ent(2),...
                lim(3,1) : lim_ent(3),...
                lim(4,1) : lim_ent(4),...
                lim(5,1) : lim_ent(5),...
                lim(6,1) : lim_ent(6),...
                ];
        idx_mem = [ lim_ent(1) : lim(1,2),...
                lim_ent(2) : lim(2,2),...
                lim_ent(3) : lim(3,2),...
                lim_ent(4) : lim(4,2),...
                lim_ent(5) : lim(5,2),...
                lim_ent(6) : lim(6,2),...
                ];

        % IDX FRECUENCIAS (MATRIZ TOTAL)
        % selecciona cada tipo de ensayo indexando la matriz total (izq/der)
        frec_tot = [ {lim(1,1):lim(1,2)} , {lim(2,1):lim(2,2)},...
                {lim(3,1):lim(3,2)}  ,  {lim(4,1):lim(4,2)},...
                {lim(5,1):lim(5,2)}  ,  {lim(6,1):lim(6,2)},...
                ];
        frec_ent = [ {lim(1,1):lim_ent(1)} , lim(2,1):lim_ent(2),...
                {lim(3,1):lim_ent(3)} , {lim(4,1):lim_ent(4)},...
                {lim(5,1):lim_ent(5)} , {lim(6,1):lim_ent(6)},...
                ];
        frec_mem = [ {lim_ent(1) :lim(1,2)},... %izq500
                {lim_ent(2) :lim(2,2)},... %izq750
                {lim_ent(3)  :lim(3,2)},... %izq1000
                {lim_ent(4) :lim(4,2)},... %der5000
                {lim_ent(5)  :lim(5,2)},... %der750
                {lim_ent(6) :lim(6,2)},...%der1000
                ];

        % % IDX STIMS
        tot_stims = [ {find(ismember(tAxis{1},0.5:0.5:4.25) )}  ;...
                {[find(tAxis{2}==0.75), find(tAxis{2}==1.5),find(tAxis{2}==2.25),find(tAxis{2}==3),...
                find( (tAxis{2} > 3.74) & (tAxis{2} < 3.76) ), find(tAxis{2}==4.5), find(tAxis{2}==5.25),...
                find( (tAxis{2} > 5.99) & (tAxis{2} <= 6.001)) ] } ;...
                {find(ismember(tAxis{3},1:1:8.5))} ;...
                {find(ismember(tAxis{1},0.5:0.5:4.25) )}  ;...
                {[find(tAxis{2}==0.75), find(tAxis{2}==1.5),find(tAxis{2}==2.25),find(tAxis{2}==3),...
                find( (tAxis{2} > 3.74) & (tAxis{2} < 3.76) ), find(tAxis{2}==4.5), find(tAxis{2}==5.25),...
                find( (tAxis{2} > 5.99) & (tAxis{2} <= 6.001)) ] } ;...
                {find(ismember(tAxis{3},1:1:8.5))} 
                ];

        ent_stims ={ [lim(1,1), find(ismember(tAxis_ent{1},0.5:0.5:1.5)) ];...
                [lim(2,1), lim(2,1) + find(ismember(tAxis_ent{2},0.75:0.75:2.25)) ];...
                [lim(3,1), lim(3,1) + find(ismember(tAxis_ent{3},1:1:3)) ];...
                [lim(4,1), lim(4,1) + find(ismember(tAxis_ent{1},0.5:0.5:1.5)) ];...
                [lim(5,1), lim(5,1) + find(ismember(tAxis_ent{2},0.75:0.75:2.25)) ];...
                [lim(6,1), lim(6,1) + find(ismember(tAxis_ent{3},1:1:3)) ];...
                };

        mem_stims ={ find(ismember(tAxis{1},1.5:.5:6));... %izq500
                (lim_ent(2)-46) + [find(ismember(tAxis{2},2.25:0.75:3)),...
                find(tAxis{2} >= 3.74 & tAxis{2} <= 3.76),...
                find(ismember(tAxis{2},4.5:0.75:5.25)),...
                find(tAxis{2} >= 5.99 & tAxis{2} <= 6.05)];... % izq750
                (lim_ent(3)-61) + find(ismember(tAxis{3},3:1:8));... %izq1000
                (lim_ent(4) - 31) + find(ismember(tAxis{1},1.5:.5:6));...  % der500
                (lim_ent(5) - 46) + [find(ismember(tAxis{5},2.25:0.75:3)),...
                find(tAxis{5} >= 3.74 & tAxis{5} <= 3.76),...
                find(ismember(tAxis{5},4.5:0.75:5.25)),...
                find(tAxis{5} >= 5.99 & tAxis{5} <= 6.05)];... % der750
                (lim_ent(6)-61) + find(ismember(tAxis{6},3:1:8)) }; % der1000


        % IDX FRECUENCIAS (MATRICES INDIVIDUALES)
        e1 = length(tAxis_ent{1}); e2 = length(tAxis_ent{2}); e3 = length(tAxis_ent{3});
        m1 = length(tAxis_mem{1}); m2 = length(tAxis_mem{2}); m3 = length(tAxis_mem{3});
        % entrainment limites [inferior, superior]
        lim_ent_only = [ lim(1,1) + e1 - 1;...
                e1 + e2 - 1;...
                e1 + e2 + e3 - 1;...
                e1*2 + e2 + e3 - 1;...
                e1*2 + e2*2 + e3 - 1;...
                e1*2 + e2*2 + e3*2 - 1;...
                ];

        frec_ent_only = { 1              : e1  ;...
                e1 + 1               : e1 + e2  ;...
                e1 + e2 + 1          : e1 + e2 + e3  ;...
                e1 + e2 + e3 + 1     : e1*2 + e2 + e3 ;...
                e1*2 + e2 + e3 + 1   : e1*2 + e2*2 + e3 ;...
                e1*2 + e2*2 + e3 + 1 : e1*2 + e2*2 + e3*2 ;...
                };

        frec_mem_only = { 1               : m1 -1;...
                m1 + 1                : m1 + m2 -1;...
                m1 + m2 + 1           : m1 + m2 + m3 -1;...
                m1 + m2 + m3 + 1      : m1*2 + m2 + m3 -1;...
                m1*2 + m2 + m3 + 1    : m1*2 + m2*2 + m3 -1;...
                m1*2 + m2*2 + m3 + 1  : m1*2 + m2*2 + m3*2 -1;...
                };    

        frec_ent_stims ={ [1,find(ismember(tAxis_ent{1},0.5:0.5:1.5)) ];...
                [lim_ent_only(1,1), lim_ent_only(1,1) + find(ismember(tAxis_ent{2},0.75:0.75:2.25)) ];...
                [lim_ent_only(2,1), lim_ent_only(2,1) + find(ismember(tAxis_ent{3},1:1:3)) ] ;...
                [lim_ent_only(3,1), lim_ent_only(3,1) + find(ismember(tAxis_ent{1},0.5:0.5:1.5)) ] ;...
                [lim_ent_only(4,1), lim_ent_only(4,1) + find(ismember(tAxis_ent{2},0.75:0.75:2.25)) ];...
                [lim_ent_only(5,1), lim_ent_only(5,1) + find(ismember(tAxis_ent{3},1:1:3))]
                };

        frec_mem_stims ={ find(ismember(tAxis_mem{1}, 2:0.5:4.25))  ;...
                [find(tAxis_mem{2} == 3), find(( tAxis_mem{2} > 3.74) & (tAxis_mem{2} < 3.76) )  ,...
                find(tAxis_mem{2} == 4.5), find(tAxis_mem{2} == 5.25)  ,...
                find(( tAxis_mem{2} > 5.98) & (tAxis_mem{2} <= 6.03) )] ;...
                find(ismember(tAxis_mem{3},4:1:8.5))
                };

end







% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% % EXAMPLE PIPELINE
% for i = 1:number_of_subjects  (brain areas)
% try
%         my_preprocessing_function(i)
%         % my_old_freqanalysis_function(i)
%         my_freqanalysis_function(i)
%         my_sourceanalysis_function(i)
%   catch
%   disp(['Something was wrong with Subject' int2str(i) '! Continuing with next in line']);
%   end
% end

% https://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/
% http://staff.www.ltu.se/~damvar/Matlab/HowToMakePrettyFiguresWithMatlab.pdf

% % FIGURE SETTINGS ------------------------------------------------------------
%     set(gca, ...
%               'Box'         , 'off'     , ...
%               'TickDir'     , 'out'     , ...
%               'TickLength'  , [.01 .01] , ...
%               'XMinorTick'  , 'off'      , ...
%               'YMinorTick'  , 'on'      , ...
%               'YGrid'       , 'on'      , ...
%               'XColor'      , [0 0 0], ...   
%               'YColor'      , [0 0 0], ...
%               'YTick'       , 0:2:maxFR, ...
%               'XTick'       , 0:1:8, ...
%               'LineWidth'   , 1         );
      
% % cambiar el color de fondo blanco/negro
% set(gca,'Color','None')
% set(gcf,'Color','None')

% ------------------------------------------------------------------------------
% % PLOT THIS

% % Stims lines 
% %  maintenance lines
% duracion = durInt();
% arrayfun(@(a)xline(a,'--','LineWidth',3),duracion*3:duracion:tAxis{3}(end)) 
% % entrainment lines
% arrayfun(@(a)xline(a,'LineWidth',3),duracion:duracion:duracion*3) 

% xlim([-0.5 tAxis{3}(end)])
% arrayfun(@(x) box(x,'off'), findobj(gcf,'Type','axes'))
% set(gcf, 'Units','normalized','Position',[ x_ y_ width hight]);

% % Common xlabel and ylabel for all subplots
% han = axes(h,'Visible','off'); % FIXME : add handle 'h' from togglefig
% han.XLabel.Visible='on'; han.YLabel.Visible='on';
% ylabel(han,'firing rate (spks/s)'); xlabel(han,'time (s)');

% % Common title
% han.Title.Visible='on';
% title(han,'yourTitle');

% ------------------------------------------------------------------------------
% para no mostrar lineas en la legend 
% legend(nameCoeffs,'Location','NorthEastOutside','AutoUpdate','off')

% ------------------------------------------------------------------------------
% Para extender valores faltantes con nan 
% a(end:valorExtendido) = missing


% ------------------------------------------------------------------------------
% ASCII code for decoration and markdown 
% ░░░░░░░░░  ▒▒▒▒▒▒▒▒▒▒  ▓▓▓▓▓▓▓▓▓
% ╔═══╦═══╗
% ║   ║   ║
% ╠═══╬═══╣
% ║   ║   ║
% ╚═══╩═══╝
% ┌──┬──┐
% │  │  │
% ├──┼──┤
% │  │  │
% └──┴──┘
% ▄▄▄▄▄▄▄▄▄▄▄ ███████████  ▀▀▀▀▀▀▀▀▀▀▀ ■■■■■■■■■■■
% ¦¦¦¦¦¦¦¦¦¦¦ ≡≡≡≡≡≡≡≡≡≡≡  ┴┴┴┴┴┴┴┴┴┴┴ ╬╬╬╬╬╬╬╬╬╬╬ 
% »»»»»»»»»»» «««««««««««  ××××××××××  °°°°°°°°°°°
% ¯¯¯¯¯¯¯¯¯¯¯
% 
%                  »»»»»»»»»»»  ver May 2022  «««««««««««
% clr
% ░░░░░░░░░░░░░░░░░░░░░░░░░░░   ~ code ~   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
% 
%      
% ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░

% ╔════════════════════════════════════════════════════════════════════════════╗
% ║ TODO:                                                                      ║ 
% 
% ╚════════════════════════════════════════════════════════════════════════════╝