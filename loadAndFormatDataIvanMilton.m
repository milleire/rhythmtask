function trials = loadAndFormatDataIvanMilton(nombreArchivo,varargin)
%
% Nueva version de loadAndFormatDataIvanMilton Jul2020.
% Usa la informacion del ensayo de otto guardada en "e.conducta" y
% en "e.anclas".
%

% Para cargar los datos con/sin LFP
incluirLFP   = getArgumentValue('incluirLFP'   ,true, varargin{:});
LFPfromDiffElectrode   = getArgumentValue('LPFfromDiffE'   ,false, varargin{:});

if incluirLFP
    graficarpWelch = 1;
else
    graficarpWelch = 0;
end

%unidad = 'I1501012D2'; % PFC
%unidad = 'I1403252A0';   % LIP Ivan

if length(nombreArchivo)==15
    nombreArchivo = nombreArchivo(1:11);
else
    nombreArchivo = nombreArchivo(1:10);
end

% Lee la estructura "n" y la estructura "e"
load(nombreArchivo);
load(strcat(nombreArchivo(1:8),'XX'));
%e.coordenadas
if incluirLFP
    % Carga el LFP y lo filtra (sampled at 1000 Hz)
    
    if LFPfromDiffElectrode % Para usar el LFP de otro electrodo
        disp('loadAndFormatDataIvanMilton.m: Usando LFP de otro electrodo')
        % Lista de archivos de ese mismo registro
%         archs=dir(['/Users/victor/Documents/mFilesData/datosIvanMilton/estructurasOtto/'...
%             nombreArchivo(1:8) '*']);
        
        % Obtiene las letras que corresponden a los electrodos
        for l=1:length(archs)
              letras(l) = archs(l).name(9);
        end
        unicas = unique(letras); % Letras unicas
        
        origLetra = nombreArchivo(9) % Original electrode
        
        if length(unicas)>1 % Si hay mas de una letra (electrodo)
            [~,locs] = ismember([origLetra 'X'],unicas);
            unicas(locs) = []; % Quita el electrodo original
        else
            disp('loadAndFormatDataIvanMilton.m: Solo hay un electrodo! Usando el mismo LFP')
            unicas = origLetra
        end
        nuevaLetra = unicas(1) % Toma la primera letra, habiendo ya quitado la original
        try
            load(strcat(nombreArchivo(1:8),[nuevaLetra 'X']))
            hayLFPs = 1;
        catch
            disp(['loadAndFormatDataIvanMilton.m: No existe archivo de LFP.'])
            hayLFPs = 0;
            c.campo = [];
        end
    else
        try
            load(strcat(nombreArchivo(1:9),'X'));
            hayLFPs = 1;
        catch
            disp(['loadAndFormatDataIvanMilton.m: No existe archivo de LFP.'])
            hayLFPs = 0;
            c.campo = [];
        end
    end

    lfpf = double(c.campo);
    if length(lfpf)>14400000
        disp('loadAndFormatDataIvanMilton.m: LFP recording too long. Trimmed at 4 Hrs.')
        disp(nombreArchivo)
        lfpf = double(c.campo(1:14400000));
    end
%    fPass = 'no filter';
    fPass = [1 110]; % Bandpass 
    [b,a] = butter(5,fPass*2/1000); 
    lfpf  = filtfilt(b,a,lfpf);
    lfpf  = int32(lfpf);
end

% Typical LFPs lengths
%   5376180
% 168304831
%  11325927
%  83674033
%   4588182
%   7350501
%   8747890 % para >450 ensayos
%   4797981
%   6295832


if (incluirLFP && graficarpWelch) && ~isempty(lfpf)
    togglefig('pwelch LFP'),clf, drawnow
    window = 2000; % window length (in samples)
    noverlap = 0;  % window overlap (in samples)
    f = 0:100;    % vector of frequencies to estimate
    [pxx, f]=pwelch(double(lfpf),window,noverlap,f,1000);

    plot(f,10*log10(pxx),'o-'), hold on
    xlabel('Frequency (Hz)'), ylabel('Power (db)')
    title([strcat(nombreArchivo(1:9),'X ') ' Filter: ' num2str(fPass)]), drawnow
end

% columnas matrizAnclas (tiempos en ms): 
% 1:inicia primero muestra
% 2:inicia periodo de memoria
% 3:tiempo del go
% 4:inicia movimiento
% 5:toca target

% Las columnas de la variable conducta son:
% (1) Tarea, 0:timing, 1:centerOut
% (2) esControl, 1:esControl, 0:noEsControl
% (3) ladoInicio, 1:izquierda, 2:derecha
% (4) duracion: 2:500ms, 3:750ms, 4:1000ms
% (5) numeroSamples
% (6) numero intervalos de memoria
% (7) esHit

nEns         = size(e.conducta,1); % numero de ensayos totales
condiciones  = double(e.conducta);
%conducta     = double(e.conducta);
matrizAnclas = double(e.anclas)/1000; % pasa de ms a s

% Tiempos de las espigas
espigas = double(n.espigas.tiempos);

eliminar = [];
% Loop over trials ----------------------------------------------------
for k = 1:nEns
    
     % Ancla es el tiempo absoluto del primero de muestra
    ancla = matrizAnclas(k,1)*1000; % we needit in ms for the LFP
    alineaEnsayo = matrizAnclas(k,1);
    
    trials(k).cEsTiming  = condiciones(k,1)==0;
    trials(k).cEsControl = condiciones(k,2);
    trials(k).cIniciaIzq = condiciones(k,3)==1;
    
    if     condiciones(k,4)==2
        trials(k).cDurIntervalo = 0.5;
    elseif condiciones(k,4)==3
        trials(k).cDurIntervalo = 0.75;
    elseif condiciones(k,4)==4
        trials(k).cDurIntervalo = 1;
    elseif condiciones(k,4)==0 % Control trials
        trials(k).cDurIntervalo = NaN;
    else
        trials(k).cDurIntervalo = NaN;
        eliminar = [eliminar k];
        disp(['loadAndFormatDataIvanMilton.m: cDurIntervalo param is not recognized.' ...
            ' Trial discarded.'])
    end
    
    trials(k).cNumeroMuestra = condiciones(k,5);
    trials(k).cNumeroMemoria = condiciones(k,6);
    trials(k).cEsHit         = condiciones(k,7);
    
    % Times of sample stimulus
    numMuestra = condiciones(k,5);
    tt = (0:numMuestra-1) * trials(k).cDurIntervalo;
    trials(k).tMuestra = tt+matrizAnclas(k,1)-alineaEnsayo;
    
    % Times of memory stimulus
    numMem = condiciones(k,6);
%    tt = (0:numMem-1)*trials(k).cDurIntervalo;
%    trials(k).tMemoria = tt+matrizAnclas(k,2)-alineaEnsayo;
    tt = (1:numMem) * trials(k).cDurIntervalo;
    if isempty(tt) || any(isnan(tt))
        trials(k).tMemoria = 0;
    else
        trials(k).tMemoria = trials(k).tMuestra(end)+tt;
    end
    
    trials(k).tGoCue  = matrizAnclas(k,3)-alineaEnsayo;
    trials(k).tIniMov = matrizAnclas(k,4)-alineaEnsayo;
    trials(k).tFinMov = matrizAnclas(k,5)-alineaEnsayo;
    trials(k).EsHit   = condiciones(k,7)==1;
    
    % lfpf is sampled at 1000 Hz so use time in ms
    lims = [-999 trials(k).tGoCue*1000];
    %trials(k)
    if isnan(lims(2))
        lims(2)=10000;
    end
    %       indices = espigas > (ancla - 1500) & espigas < ancla + 10000;
    indices = espigas > (ancla - 2000) & espigas < (ancla + lims(2)+1000);
    %   indices = espigas > ancla - 1000 & espigas < anclaFin;
    trials(k).spikeTimes = double(espigas(indices) - ancla) / 1000;
    
    if incluirLFP && ~any(ismember(eliminar,k))
%    if incluirLFP 
        
        if isfield(c,'campoTimeStamps') % For the files with pauses
            tIni = (ancla+lims(1))/1000;
            tFin = (ancla+lims(2))/1000;
            [~, idxIni] = min(abs(c.campoTimeStamps - tIni));
            [~, idxFin] = min(abs(c.campoTimeStamps - tFin));
        else
            idxIni = round(ancla+lims(1));
            idxFin = round(ancla+lims(2));
        end
        
        if idxFin < length(lfpf)
            trials(k).LFP   = single( lfpf( idxIni : idxFin ));
        else
            trials(k).LFP   = [];
        end
        
        trials(k).LFPt0 = -lims(1) + 1; % Position of the time=0
        
    end
    
    % For the task trials, get the side of the screen that was touched
%    if trials(k).cEsTiming && ~trials(k).cEsControl  
%        trials(k).cIniciaIzq
    
%    trials(k).cDecisionIsLeft = 
    
%    trials(k)
%    keyboard
end

trials(eliminar)=[];












