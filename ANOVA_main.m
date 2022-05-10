%              ANOVA para clasificacion neuronas de hipocampo
%               »»»»»»»»»»»  ver May 2022  «««««««««««

clr
% ░░░░░░░░░░░░░░░░░░░░░░░░░░░   ~ code ~    ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
%
% orden de analisis : 
% 1) Se crearon neuronas modeladas
%   * Especificaciones * 
%   - las de preferencia lado izq, inician los trials en lado izquierdo
%   - N neurons per matrix
%   - Classes :  
%       - ramping
%       - max activity on 2nd interval
%       - left side preference 
%       - absolute time on 750 ms
%
% 2) Se realizo analisis anova con neuronas reales y modeladas
%       - Factores del analisis
%           fr       : meanfr de cada ventana de análisis
%           sideStim : lado del estimulo presentado (categorico)
%           sideInit : lado del estimulo de inicio del trial (categorico)
%           dur      : duración del estimulo en el trial (continuo)
%           interv   : intervalo en el que se encuentra la ventana (continuo)
%       - Interacciones relevantes
%           durs * SideStim
%           durs * Interval
%           SideStim * SideInit
%           SideStim * Interval
% 
% 3) Se plotearon los resultados de las diferentes subpoblaciones.
%   - se clasifico a cada neurona dependiendo de su pval significativo para
%       cada factor
% 
% NOTA:
% - El primer analisis se hizo con la FR de fr_hipp_202204.mat.  El
% siguiente se hizo con fr_hipp_202205.mat  Se cambió la ventana de fr a
% 0.3ms
% 
% ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░

% ╔════════════════════════════════════════════════════════════════════════════╗
% ║ TODO:                                                                      ║ 
% - ifWraping :eje de tiempo comun. linspace. variar ancho de bines
% - modificar la timeConstant, se ve muy ruidosa la fr 
% ╚════════════════════════════════════════════════════════════════════════════╝



% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
%                             NEURONAS MODELADAS
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
% 'ramp_matrix',  'max2_matrix', 'side_matrix', 'abs75_matrix' 
%  Cada una tiene 3x1 cell  ({1}:500ms ; {2}:750ms ; {3}:1000ms)
%  En cada celda tienen 150 neuronas con 200 trials cada una 
area = 1;
run('main.m')

% variables
N = 144; % neuronas x clase 
T = 200; % trials x neurona
minval = 0.01;
maxval = 1.5;
subpp  = repmat(1:3, [1,2]);

% matrices que almacenan las 3 duraciones. Cada celda contiene matrices : 
ramp_matrix = cell(N,6); % (N, size(tAxis{ii},2))
max2_matrix = cell(N,6); % (N, size(tAxis{ii},2))

% para side_matrix N*2 en adelante son preferencia lado derecho
side_matrix = cell(N,6);  % (N*2, size(tAxis{ii},2))  
abs75_matrix = cell(N,6); % (N, size(tAxis{ii},2))

% condicionales
graficar = [ 0 0 0 0 ];  % ramping, max2nd, leftPref, abs750
stats    = [ 0 0 0 0 ];
saveVars = 1;


%...............................................................................
% RAMPING NEURONS --------------------------------------------------------------
k = 0.1;
init = [ 5 5 5 5 5 5 ];

% construye matriz de neuronas por cada duracion
for ii = 1:6 % loops durs (.5 .75 1 .5 .75 1)
    x = tAxis{ii};
    for jj = 1:N % loops number of neurons 
        cumulativeTT = nan(T,size(x,2));
        for tt = 1:T
            % suma noise más fr 
            noise = -1 + minval + rand(1,size(x,2)) .* (maxval - minval);
            a = init(ii) * exp(k .* x );
            cumulativeTT(tt,1:size(a,2)) =  a + noise;
        end
       ramp_matrix{jj,ii} = cumulativeTT; 
    end
end



if graficar(1)
    for ii = 1:6
        x = tAxis{ii};
        
        togglefig('ramping cells')
        meanfr = mean(cell2mat(ramp_matrix{1,ii}),1);
        stdfr  = std(cell2mat(ramp_matrix{1,ii}) );
        subplot(3,1,subpp(ii))
        plot(x,meanfr, 'LineWidth' , 3,'color', 'k'), hold on
        patch([x, fliplr(x)],[meanfr-stdfr, fliplr(meanfr+stdfr)], 'g',...
                    'FaceAlpha',0.3, 'EdgeColor','none'); 
                
        arrayfun(@(a)xline(a,'--','LineWidth',4),durInt(ii)*3:durInt(ii):x(end)) 
        arrayfun(@(a)xline(a,'LineWidth',3),durInt(ii):durInt(ii):durInt(ii)*3) 
        xlim([0 tAxis{3}(end)])
        box off
        set(gcf, 'Units','normalized','Position',[0 0.5 0.5 0.4]);        
    end
end



if stats(1)
    % ttest
    for ii = 1:6
        x = tAxis{ii};
        % idx 1er intervalo
        idx1 = find(x== durInt(ii))   : find(x== durInt(ii)*2) ;
        % idx 2do intervalo
        idx2 = find(x== durInt(ii)*7) : find(x==durInt(ii)*8)  ; 
        fr = cell2mat(ramp_matrix{1,ii});
        x_data = mean(fr(:,idx1),2);
        y_data = mean(fr(:,idx2),2);

        [h, p] = ttest(x_data,y_data);
        if any( p < 0.001) 
            disp(['significativo para ramping activity : dur ',...
                num2str(durInt(ii))])
        end
    end
end



%...............................................................................
% MAX ACTIVITY 2nd DURATION ----------------------------------------------------
k    = [3 1.5 1 3 1.5 1];
init = [60 40 60 60 40 60];
% construye matriz de neuronas cada duracion
for ii = 1:6 
    x = tAxis{ii};
    
    % calcula la exponencial CRECIENTE antes de 2do interval
    a = exp( x(1:find(x == durInt(ii)*2)) );
    
    % calcula la exponencial DECRECIENTE después de 2do interval
    b = init(ii) .* exp(- k(ii) .* x(find(x == durInt(ii)*2) + 1:end) );

    for jj = 1:N
        cumulativeTT = nan(T,size(x,2));
        for tt = 1:T
            % suma noise más fr 
            noise = -1 + minval + rand(1,size(x,2)) .* (maxval - minval);
            cumulativeTT(tt,:) = [ a , b ] + noise;
        end
        max2_matrix{jj,ii} = cumulativeTT;
    end
    
end

if graficar(2)
    togglefig('maximo en segundo intervalo')
    for ii = 1:6
        x = tAxis{ii};
        meanfr = mean(cell2mat(max2_matrix{1,ii}),1);
        stdfr  = std(cell2mat(max2_matrix{1,ii}));
        
        subplot(3,1,subpp(ii))
        plot(x,meanfr, 'LineWidth' , 3 ,'color', 'k'), hold on
        patch([x, fliplr(x)],[meanfr-stdfr, fliplr(meanfr+stdfr)], 'g',...
                    'FaceAlpha',0.3, 'EdgeColor','none'); 
        
        arrayfun(@(a)xline(a,'--','LineWidth',4),durInt(ii)*3:durInt(ii):x(end)) 
        arrayfun(@(a)xline(a,'LineWidth',3),durInt(ii):durInt(ii):durInt(ii)*3) 
        xlim([0 tAxis{3}(end)])
        box off
    end
    set(gcf, 'Units','normalized','Position',[0 0 0.5 0.4]);
end

if stats(2)
    % ttest
    
    for ii = 1:3
        x = tAxis{ii};
        % idx actividad maxima en 2do intervalo
        if ismember(ii, [1,3])
            idx1 = find(x==durInt(ii)*1.5) : find(x== durInt(ii)*2.5) ;
        else 
            idx1 = find(x==1.12) : find(x==1.87) ;
        end
        
        % idx actividad basal
        idx2 = find(x== durInt(ii)*3) : find(x==durInt(ii)*8);  
        fr = cell2mat(max2_matrix{1,ii});
        x_data = mean(fr(:,idx1),2);
        y_data = mean(fr(:,idx2),2);
        
        [h, p] = ttest(x_data,y_data);
        if any( p < 0.001) 
            disp(['significativo para max activity 2do interval : dur ',...
                num2str(durInt(ii))])
        end
    end
end








%...............................................................................
% PREFIERE LADO IZQUIERDO ------------------------------------------------------
amplitude      = [5, 6, 5.5, 5, 6, 5.5]  ;
period         = [6.4, 4.2, 3.1, 6.4, 4.2, 3.1]; 
minval = 0.01;
maxval = 3;

% construye matriz de neuronas cada duracion
for ii = 1:6
    x = tAxis{ii};
    for jj = 1 : N
        cumulativeTT = nan(T,size(x,2));
        for tt = 1:T        
            if ismember(ii,1:3)
                % preferencia lado izquierdo
                phase_shift    = [-0.29, -0.4, -0.5, -0.29, -0.4, -0.5];  
                noise   = 1+ minval + rand(1,size(x,2)) .* (maxval - minval);
                cumulativeTT(tt,:) = amplitude(ii) + sin(period(ii) *  ...
                    (x + phase_shift(ii)) ) + noise;
            elseif ismember(ii,4:6)
                % preferencia lado derecho
                phase_shift     = [-0.8, -1.2, -1.5, -0.8, -1.2, -1.5];  
                vertical_shift  = -3;
                noise   = 1 + minval + rand(1,size(x,2)) .* (maxval - minval);
                cumulativeTT(tt,:) = amplitude(ii) + sin( period(ii) * ...
                    (x + phase_shift(ii)) ) + vertical_shift + noise;
            end
        end
        side_matrix{jj,ii} = cumulativeTT;
    end    
end


if graficar(3)
    togglefig('prefieren lado izquierdo')
    for ii = 1:6
        x = tAxis{ii};
        
        meanfr = mean(cell2mat(side_matrix{1,ii}),1);
        stdfr  = std(cell2mat(side_matrix{1,ii}));
        
        subplot(3,1,subpp(ii))
        % lado izquierdo
        plot(x,meanfr, 'LineWidth' , 5 ,'color', c.linea{ii}), hold on
        patch([x, fliplr(x)],[meanfr - stdfr, fliplr(meanfr + stdfr)],...
            c.linea{ii}, 'FaceAlpha',0.2, 'EdgeColor','none'); 
        
        arrayfun(@(a)xline(a,'--','LineWidth',4),durInt(ii)*3:durInt(ii):x(end))
        arrayfun(@(a)xline(a,'LineWidth',3),durInt(ii):durInt(ii):durInt(ii)*3) 
        xlim([0 tAxis{3}(end)])
        box off
        set(gcf, 'Units','normalized','Position',[0.5 0.5 0.5 0.4]);
    end
end



if stats(3)
    for ii = 1:3
        % ttest
        x = tAxis{ii};
        idx = 1 : find(x==2.5);

        % lado izquierdo
        data = cell2mat(side_matrix{1,ii});
        x_data = data(:,idx );

        % lado derecho
        data = cell2mat(side_matrix{1,ii+3});
        y_data = data(:,idx);

        [h, p] = ttest(x_data,y_data);
        if sum( p < 0.001) == size(x_data,2)
            disp(['significativo para Left vs Right, para las durs : ',...
                num2str(durInt(ii))])
        end
    end
end



%...............................................................................
% % tiempo absoluto en 750ms (sin importar frecuencia) -------------------------
minval = 0.1;
maxval = 3;

k      = [ 2  2  2  2  2  2 ];
initA  = [ 1  1  1  1  1  1 ];
initB  = [ 20 20 20 20 20 20 ];

for ii = 1:6
    x = tAxis{ii};

    for jj = 1:N
        cumulativeTT = nan(T,size(x,2));
        for tt = 1:T    
            a = initA(ii) * exp( k(ii) .*  x(1:find(x == 0.75)) );
            b = initB(ii) .* exp(-k(ii) .* x(find(x == 0.75) + 1 : end) );
            temp = [ a , b ];
            noise = minval + rand(1,size(temp,2)) .* (maxval - minval);
            cumulativeTT(tt,1:size(temp,2)) = [ a , b ] + noise;
        end
        abs75_matrix{jj,ii} = cumulativeTT;
    end
    
end

if graficar(4)
    togglefig('Absoluto en 750')
    for ii = 1:6
        x = tAxis{ii};
        meanfr = mean(cell2mat(abs75_matrix{1,ii}),1);
        stdfr  = std(cell2mat(abs75_matrix{1,ii}));

        subplot(3,1,subpp(ii))
        plot(x,meanfr, 'LineWidth' , 3 ,'color', 'k'), hold on
        patch([x, fliplr(x)],[meanfr-stdfr, fliplr(meanfr+stdfr)], 'g',...
                    'FaceAlpha',0.3, 'EdgeColor','none'); 
        
        arrayfun(@(a)xline(a,'--','LineWidth',4),durInt(ii)*3:durInt(ii):x(end)) 
        arrayfun(@(a)xline(a,'LineWidth',3),durInt(ii):durInt(ii):durInt(ii)*3) 
        xlim([0 tAxis{3}(end)])
        box off
        
    end
    set(gcf, 'Units','normalized','Position',[0.5 0 0.5 0.4]);
end


if stats(4)
    for ii = 1:6
        
        % ttest
        x = tAxis{ii};
        idx1 = find(x == 0 ) : find(x == 1.5 ) ;
        idx2 = find(x == durInt(ii)*3) : find(x == durInt(ii)*8) ;

        % lado izquierdo
        fr = cell2mat(abs75_matrix{1,ii});
        x_data = mean(fr(:,idx1), 2);
        y_data = mean(fr(:,idx2), 2);

        [h, p] = ttest(x_data,y_data);
        
        if any(p < 0.001) 
            disp(['significativo para max activity on 750ms : ',...
                num2str(durInt(ii))])
        end
        
    end
end


%...............................................................................
% ------------------------------------------------------------------------------
FR_dummy = [ramp_matrix; max2_matrix; side_matrix; abs75_matrix(1:N-2,:) ];

if saveVars
    % saves all neurons 
    save(fullfile(mainpath,'dummyfr_hipp_202205.mat'), 'FR_dummy', '-v7.3')  
end

disp('finalizado FR dummy')














%%
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
%                              ANOVA ANALYSIS 
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
clr
area = 1;
run('main.m')

% conditionals
saveVars = 1;
inicia = 1;

% load data 
load('dummyfr_hipp_202205.mat')
load('fr_hipp_202205.mat')

% anova term and their interactions
sources = {'durs','sideStim','sideInit','Interv','durs*sideStim','durs*Interv',...
                'sideStim*sideInit','sideStim*Interv','durs*sideStim*Interv'};
            
            
% Preallocation and var definition 

% struct neuron para guardar datos de resultados de c/neurona
neuron       = struct('pval',[],'stats',[],'tbl',[]);
neuron.pval  = zeros(size(FR,1),size(sources,2));
neuron_dummy = neuron;

% struct anova para guardar factores del analisis de c/neurona
anovamat       = struct('fr',[],'sideStim',[],'sideInit',[],'dur',[],...
                'interv',[],'epoch',[]);
anovamat_dummy = anovamat;

numBines = 100;
tAxisWid = [ {linspace(0,tAxis{1}(end),numBines )},...
             {linspace(0,tAxis{2}(end),numBines )},...
             {linspace(0,tAxis{3}(end),numBines )},...
             {linspace(0,tAxis{1}(end),numBines )},...
             {linspace(0,tAxis{2}(end),numBines )},...
             {linspace(0,tAxis{3}(end),numBines )} ];

idxStims = [ {durInt(1):durInt(1):durInt(1) * 9},...
             {durInt(2):durInt(2):durInt(2) * 9},...
             {durInt(3):durInt(3):durInt(3) * 9},...
             {durInt(1):durInt(1):durInt(1) * 9},...
             {durInt(2):durInt(2):durInt(2) * 9},...
             {durInt(3):durInt(3):durInt(3) * 9} ];
             
esLado   = [ {repmat([IniciaLado(1) , ~IniciaLado(1)],[1,5])},...
             {repmat([IniciaLado(2) , ~IniciaLado(2)],[1,5])},...
             {repmat([IniciaLado(3) , ~IniciaLado(3)],[1,5])},...
             {repmat([IniciaLado(1) , ~IniciaLado(1)],[1,5])},...
             {repmat([IniciaLado(2) , ~IniciaLado(2)],[1,5])},...
             {repmat([IniciaLado(3) , ~IniciaLado(3)],[1,5])} ];
   
         

for nn = inicia:size(FR,1)  % nn = neuronas
    savetrials = [];
    cs  = [ {1},{1},{1},{1},{1},{1} ];
    SZ  = 1;
    SZd = 1;    
    
    
    for tt = 1:length(tAxisWid{1})-1  % tt = window        
        
        for jj = 1:length(IniciaLado) % jj = trials
            binInicia  = tAxisWid{jj}(tt);
            binFin     = tAxisWid{jj}(tt+1);
            
            % selection of window
            idx = ( tAxis{jj} >= binInicia ) & ( tAxis{jj} <= binFin );
            
            temp    = mean(FR{nn,jj}(:,idx), 2, 'omitnan' ); 
            sz      = size(temp,1);
            
            tempD   = mean(FR_dummy{nn,jj}(:,idx), 2, 'omitnan' ); 
            szd     = size(tempD,1); 
            
            
            % saves factors for analysis ..................
            
            % FIRING RATE
            anovamat.fr( SZ:sz + SZ -1, 1)                = temp;
            anovamat_dummy.fr( SZd:szd + SZd -1, 1)       = tempD;
            
            % SIDEINIT
            anovamat.sideInit( SZ:sz + SZ -1, 1)          =...
                logical( repmat(IniciaLado(jj),[sz,1]) );
            anovamat_dummy.sideInit( SZd:szd + SZd -1, 1) =...
                logical( repmat(IniciaLado(jj),[szd,1]) );
            
            % DUR
            anovamat.dur( SZ:sz + SZ -1, 1)              =...
                repmat(durInt(jj),[sz,1]);
            anovamat_dummy.dur( SZd:szd + SZd -1, 1)     =...
                repmat(durInt(jj),[szd,1]);
            
            % INTERV
            timePoint = ( binFin - binInicia ) / 2;
            anovamat.interv( SZ:sz + SZ -1, 1)           =...
                repmat(timePoint,[sz,1]);
            anovamat_dummy.interv( SZd:szd + SZd -1, 1)  =...
                repmat(timePoint,[szd,1]);
            
            % SIDESTIM
            anovamat.sideStim(SZ:sz + SZ -1, 1)          = ...
                    repmat(esLado{jj}(cs{jj}),[sz,1]);
                    
            anovamat_dummy.sideStim(SZd:szd + SZd -1, 1) = ...
                    repmat(esLado{jj}(cs{jj}),[szd,1]);
                
            
            if abs(binFin - idxStims{jj}(cs{jj}) ) < 0.08
                cs{jj} = cs{jj} + 1;
            end
            
            
            % EPOCH
            if cs{jj} <= 3
                anovamat.epoch(SZ:sz + SZ -1, 1)          = ones(sz,1);
                anovamat_dummy.epoch(SZd:szd + SZd -1, 1) = ones(szd,1);
            else
                anovamat.epoch(SZ:sz + SZ -1, 1)          = ones(sz,1)*2;
                anovamat_dummy.epoch(SZd:szd + SZd -1, 1) = ones(szd,1)*2;
            end
                
                         
            SZ  = SZ + sz;
            SZd = SZd + szd;
            
            
        end
        
        
        terminos = [ 1 0 0 0 0 ; 1 1 0 0 0 ; 1 1 1 0 0 ; 1 1 1 1 0 ; 1 1 1 1 1];

        % ANOVA reales ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        [~,~,stats,~] = anovan(anovamat.fr,...
            {anovamat.dur,anovamat.sideStim,anovamat.sideInit, anovamat.interv,...
            anovamat.epoch},...
            'alpha', 0.05,...
            'Varnames',{'durs','sideStim', 'sideInit', 'Interv', 'Epoch'},...
            'Continuous', [1, 4, 5],...
            'model', terminos,...
            'Display','off'); 

%         % pval significativo -> categorias 
%         if any( p < 0.05 )
%             idx = find(p < 0.05);
%             names = tbl(idx+1,1);
% 
%             for cc = 1:size(sources,2)
%                 neuron.pval(nn,cc) = any(strcmp(names,sources{cc}));
%             end
% 
%         end

        % saves results
        neuron.stats{nn,tt}         = stats.coeffs;
%         neuron.tbl{nn}           = tbl;


        % ANOVA dummy ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
        [~,~,stats,~] = anovan(anovamat_dummy.fr,...
            {anovamat_dummy.dur,anovamat_dummy.sideStim,anovamat_dummy.sideInit,...
            anovamat_dummy.interv, anovamat_dummy.epoch},...
            'alpha', 0.05,...
            'Varnames',{'durs','sideStim', 'sideInit', 'Interv', 'Epoch'},...
            'Continuous', [1, 4, 5],...
            'model', terminos,...
            'Display','off'); 

%         % pval significativo -> categorias 
%         if any( p < 0.05 )
%             idx = find(p < 0.05);
%             names = tbl(idx+1,1);
% 
%             for cc = 1:size(sources,2)
%                 neuron_dummy.pval(nn,cc) = any(strcmp(names,sources{cc}));     
%             end
%         end


        % saves results
        neuron_dummy.stats{nn,tt}         = stats.coeffs;
%         neuron_dummy.tbl{nn}           = tbl;


    end
    
    
    
end

% Saves analysis data struct  
if saveVars
    save(fullfile(mainpath,['anova_',brainAreas{area},'.mat']), 'sources',...
        'neuron','neuron_dummy')
end

disp('finished ANOVA')


%%
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
%                               PLOT RESULTS 
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 

% clr
clc
area = 1;
run('main.m')

% % load data
% load('fr_hipp_202205.mat')        % FR de neuronas reales
% load('dummyfr_hipp_202205.mat')   % FR neuronas dummy
% load('anova_hipp.mat')            % resultados ANOVA FR real & dummy

% ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
% Files contain : 
% ■ 'fr_hipp_202204.mat' : 'FR' neuronas reales. cell(Neurons x 6)
% ■ 'dummyfr_hipp_202204.mat' : 'FR_dummy' neuronas modelo. cell(Neurons x 6), 
%           y matrices individuales : 'ramp_matrix','max2_matrix','side_matrix',
%          'abs75_matrix'
% ■ 'anova_hipp.mat': estructuras neuron, y neuron_dummy. 
%           Con fields 
%               - pval, (Neurons x Num de sources) indice logico de sources con
%                  pval < 0.05
%               - stats, resultado 'stats' de la funcion anovan de c/neurona
%               - tbl, resultado 'tbl' de la funcion anovan de c/neurona
% ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

neur = 2;
tAxisWid = [ {linspace(0,tAxis{1}(end),numBines )},...
             {linspace(0,tAxis{2}(end),numBines )},...
             {linspace(0,tAxis{3}(end),numBines )},...
             {linspace(0,tAxis{1}(end),numBines )},...
             {linspace(0,tAxis{2}(end),numBines )},...
             {linspace(0,tAxis{3}(end),numBines )} ];
subpp = [ 1 2 3 1 2 3 ];

% grafica la fr
togglefig('firing rate')
for jj = 1:6
    subplot(3,1,subpp(jj))
    plot(tAxis{jj}, mean(FR_dummy{neur,jj},1), 'Color', c.redblue{jj}, 'LineWidth', 5)
    hold on 
    
    % Stims lines 
    %  maintenance lines
    duracion = durInt(jj);
    arrayfun(@(a)xline(a,'--','LineWidth',3),duracion*3:...
        duracion:tAxis{3}(end)) 
    % entrainment lines
    arrayfun(@(a)xline(a,'LineWidth',3),duracion:duracion:...
        duracion*3) 
    xlim([-0.5 tAxis{3}(end)])
end


togglefig('coefficients')
for ii = 1 : 99
    scatter(ii*ones(1,size(neuron_dummy.stats{neur,ii},1)),...
        neuron_dummy.stats{neur,ii})
    hold on
end




























%%
% % BACKUP
% % % Condicionales 
% % graficar    = 1;
% % 
% % % Definicion de variables
% % N           = 4; % num de neuronas para visualizar x categoria
% % numsubplots = 3*N;
% % subpp       = [ reshape(1:numsubplots, [N,3])' ; reshape(1:numsubplots, [N,3])' ];
% % structNames = {neuron, neuron_dummy}; % database anova results
% % getFR       = {FR , FR_dummy}; %database
% % tags        = {'neuronas reales', 'neuronas modelo'}; % database names
% % 
% % % Posibles combiaciones de factores anova 
% % combinaciones = (dec2bin(0:2^9-1)' - '0')';
% % saveCat       = zeros(size(combinaciones,1),2);
% % 
% % a = 0;
% % for db = 2 : 2 %numel(structNames)  % 1:reales  2:dummy  
% %     S             = structNames{db};
% %     firingRateMat = getFR{db};
% %     
% %     % Define clasificaciones a partir de las sources significativas en pval
% %     idx = find(sum(S.pval,1)>1);  
% %     factorNames = sources(idx);
% %     % 9 factores -> numel(idx) factores
% %     idxPval = unique( combinaciones(:,idx),'rows' );  
% %     
% %     for ss = 1 : size(idxPval,1)
% %         
% %         % Selecciona neuronas con combinaciones de factores significativos
% %         A = repmat(idxPval(ss,:), size(S.pval,1),1);
% %         temp = arrayfun(@ismember, A, S.pval(:,idx));
% %         idxFR = all(temp,2);
% %         
% %         if ~sum(idxFR)==0
% %             fr = firingRateMat(idxFR,:);
% %             
% %             catName = '';
% %             for ee = 1:size(idxPval(ss,:),2)
% %                 if idxPval(ss,ee) == 1
% %                     catName = [catName,'/ ', factorNames{ee},  ' /' ];
% %                 end
% %             end
% %             
% %             if isempty( catName )
% %                 continue
% %             end
% %             
% %             disp('°°°°°°°°°°°°°°°°°°°°°°')
% %             disp(catName)
% %             disp(['numero de neuronas :', num2str(size(fr,1)) ])
% %             a = a + size(fr,1);
% %             if graficar
% %                 
% %                 
% %                 h = togglefig(catName);
% %                 
% %                 for pot = 1:N
% %                     fr_ = fr(randi(size(fr,1)),:);            
% %                     for d = 1:numel(durInt)
% %                         
% %                         subplot(3,N,subpp(d,pot) );
% %                         plot(tAxis{d},fr_{d},...
% %                             'LineWidth', 3,...
% %                             'color', c.lineaLR{d})
% %                         hold on
% %                         
% %                         % Stims lines 
% %                         %  maintenance lines
% %                         duracion = durInt(d);
% %                         arrayfun(@(a)xline(a,'--','LineWidth',3),duracion*3:...
% %                             duracion:tAxis{3}(end)) 
% %                         % entrainment lines
% %                         arrayfun(@(a)xline(a,'LineWidth',3),duracion:duracion:...
% %                             duracion*3) 
% %                         xlim([-0.5 tAxis{3}(end)])
% %                     end
% %                 end
% % 
% %                 han = axes(h,'Visible','off');
% %                 han.XLabel.Visible='on'; han.YLabel.Visible='on';
% %                 ylabel(han,'firing rate (spks/s)'); xlabel(han,'time (s)');
% %                 han.Title.Visible='on'; title(han,[tags{db}, ' -> Categoria: ',...
% %                     catName])
% %                 arrayfun(@(x) box(x,'off'), findobj(gcf,'Type','axes'))
% %                 set(gcf, 'Units','normalized','Position',[ 0 0 .8 .85]);
% %             end
% %             
% %         end
% %         
% %     end
% % end
% % 
% % 
% % % TODO: guardar las fr de cada categoria en pdfs separados 
% % % hacer una gráfica de pastel con las clasificaciones.. o un venn diagram





%% 
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 
%                                  REFERENCES
% ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ 



% ejemplos y documentacion
% https://www.mathworks.com/help/stats/n-way-anova.html
% http://compneurosci.com/wiki/images/e/e6/Repeated_ANOVA_MATLAB_v2.pdf

% multiple comparisons
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6193594/
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7720730/

% interpret multcompare
% https://www.cmaj.ca/content/166/1/65.short

% regresores anovan
% https://www.mathworks.com/help/stats/f-statistic-and-t-statistic.html
% https://www.mathworks.com/help/stats/linearmodel.anova.html
% https://www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html


