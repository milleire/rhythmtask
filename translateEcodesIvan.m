function translation = translateEcodesIvan(codesToTranslate)
%
% Examples:
% ecode = translateEcodes('saccadeini')
% ecode = translateEcodes(7120)
% ecode = translateEcodes('spike11')   % code of a spike in channel 1, unit 1.
% ARREGLAR UN BUG: translateEcode(7)

% Translation table. (TS: timestamp)
% lookUpTable = [...
% {'fixationCentOut'}       {'18'};
% {'delayGoMov'}            {'19'};
% {'interTrialCentOut'}     {'20'};
% {'hitCentOut'}            {'21'};
% {'error'}                 {'22'};
% {'touchTarHit'}           {'23'};
% {'goMovTouchTar'}         {'24'};
% {'memRight'}              {'25'};
% {'memLeft'}               {'26'};
% {'samRight'}              {'27'};
% {'samLeft'}               {'28'};
% {'fixation'}              {'29'};
% {'interTrial'}            {'30'};
% {'spike'}               {'9000'}];

% Para el re-entramiento de Ivan 2018
lookUpTable = [...
{'fixationCentOut'}       {'18'};
{'delayGoMov'}            {'19'};
{'interTrialCentOut'}     {'20'};
{'hitCentOut'}            {'21'};
{'error'}                 {'22'};
{'touchTarHit'}           {'23'};
{'goMovTouchTar'}         {'24'};
{'memRight'}              {'25'};
{'memLeft'}               {'26'};
{'samRight'}              {'27'};
{'samLeft'}               {'28'};
{'fixation'}              {'29'};
{'interTrial'}            {'30'};
{'spike'}               {'9000'}];

% En el center out: Del 24 al 29 son las posiciones. El 24 es 0 grados, 25 es 60

% Allocate space for the translated output
if iscell(codesToTranslate)
   translation = zeros(size(codesToTranslate));
elseif ischar(codesToTranslate) %If input is a single string
   code                = codesToTranslate;
   codesToTranslate    = [];
   codesToTranslate{1} = code;
else
   translation = cell(size(codesToTranslate));
end

% Loop through each input
for k =  1 : length(codesToTranslate)   
   if iscell(codesToTranslate)    % If input is a string translate it to a number
%      codeToTranslate = lower(codesToTranslate{k});
       codeToTranslate = codesToTranslate{k};
      if length(codeToTranslate)>=5 && strcmp(codeToTranslate(1:5),'spike') % If string is a spike name
         translation(k) = 9000+(str2double(codeToTranslate(6:end)));
      else % String must be an ecode, not a spike name
         index = strmatch(codeToTranslate,lookUpTable(:,1),'exact');
         if isempty(index)
            disp([mfilename '.m :: Ecode ' codesToTranslate{k}  ' not found. Returning -1.'])
            translation(k) = -1;
         else
            translation(k) = str2double(lookUpTable{index,2});
         end
      end
   else    % Else, input must be a number so translate it to a string
      if codesToTranslate(k) >=9000 % The number is a spike name
         translation{k} = ['spike' num2str(codesToTranslate(k)-9000)]; % Create spike name
      else  % The number is an ecode
         index = strmatch(num2str(codesToTranslate(k)),lookUpTable(:,2),'exact');
         if isempty(index)
             %disp([mfilename '.m :: Ecode ' num2str(codesToTranslate(k))  ' not found. Returing number as string.' ])
             translation{k} = ['   ' num2str(codesToTranslate(k))];
         else
             translation{k} = lookUpTable{index,1};
         end
      end
      if length(translation) == 1  % If input is a single number, return a string, not a cell containing the string.
         translation = translation{1};
      end
   end
end

