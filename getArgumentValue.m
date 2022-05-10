function value = getArgumentValue(argument,defaultValue,varargin)
% getArgumentValue() Is used within a function to get the value following an
%                    argument in an 'argument',[value] pair within VARARGIN.
%
%   Loops through VARAGIN to find ARGUMENT and get its value. If ARGUMENT is not 
%   found VALUE equals DEFAULTVALUE.
%
%   Use:
%   value = getArgumentValue('argument',defaultvalue,varargin{:})
%                     Use defaultvalue if 'argument' is not found.
%
%   value = getArgumentValue('argument',defaultvalue,varargin{:},'warningoff')
%                     Use defalut value, but don't issue a warning. 
%
%   []    = getArgumentValue('argument',[],varargin)
%                     No default value is necessary.
%
%   []    = getArgumentValue('argument','msg:warning msg',varargin{:})
%                     No default value is necessary, but display a msg.
%
%   []    = getArgumentValue('argument','errorIfEmpty',varargin{:})
%                     The argument-value pair must be provided.
%
%   Examples:
%      windowLen = getArgumentValue('window'   , 0.5       , varargin{:})
%      align     = getArgumentValue('align'    , 'dots_on' , varargin{:})
%      quickplot = getArgumentValue('quickplot', 'yes'     , varargin{:})
%
% vdl 2007 - 2020

% Loop through VARARGIN searching for the ARGUMENT string.
found=0; k=1;
while ~found && k<=length(varargin)
   if ischar(varargin{k})
      if strmatch( lower(argument) , lower(varargin{k}), 'exact' )
         found = k;
      end
   end
   k=k+1;
end

if found
   value = varargin{found+1}; % Return the value, following the ARGUMENT string.
else
   value = defaultValue;         % If ARGUMENT not found return the defaultValue. 

   callF = dbstack;              % Get the names of the calling functions. 
   indx  = (length(callF)>1)+1;  % Index to the name of the calling function.

   % Issue an error if the 'errorIfEmpty' string is found in the DEFAULTVALUE.
   if ischar(value) && (strcmpi(value,'errorIfEmpty'))
      error(['Input argument ''' argument ''' for ' callF(indx).name '() not found.'])

   % Return an empty value but issue a warning if the string 'msg:' is found in the DEFAULTVLUE.
   elseif ischar(value) && length(value)>3 && ~isempty(strmatch('msg:',value))
      value=[];
      %disp([callF(indx).name '() ' defaultValue(5:end)])

   % Display a reminder message if the string 'warningOff' is not found as the last element of varargin.
%   elseif ~isempty(value) && (isempty(varargin)  || isempty(strmatch('warningoff',lower(varargin{end}))))
   elseif ~isempty(value) && (isempty(varargin)  )
      %disp(['Input argument ''' argument ''' for ' callF(indx).name '() not found. Using default value: ' num2str(value)])
   end
end
