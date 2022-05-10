function trials = selTrialsIvanMilton(trials,varargin)
%
%  Examples:
%  selected = selectTrials(trials(1:100),'corr_tarangle',[0 180],'selectneuron','spike21') % Spikes from channel 2, unit 1.
%  selected = selectTrials(trials(1:100),'align','dots_on','corr_tarangle',[270 90])       % Note that [LowerLim > UpperLim].
%
%  Note: [Lower Upper] limits define closed intervals, i.e. they include the limit values.
%        if     [Lower > Upper], the selected values are [Lower:+Inf & -Inf:Upper],
%        elseif [Lower < Upper], the selected values are [Lower : Upper]
%  'prefDir',[dir] Uses the neuron's preferred direction (0 or 180) to set the trial.direction field to 0 (nonPref) or 1 (prefDir).


n_trials     = length(trials);                        % Number of trials.
RemoveTrials = zeros(n_trials, 1);                    % Trials that will be removed after the selection.
fields       = [fieldnames(trials);'dotsduration'];   % Obtain field names.
alignEvent   = [];                                    % The default is do not aligment.


% For selecting those trials following a transition from one task to another. 
% transTrials = getArgumentValue('transitionTrials',[],varargin{:}); % Proportion of trials after a transition that get selected
% if ~isempty(transTrials)
%    selectID = getArgumentValue('transitionTrials','errorifempty',varargin{:});
%    trialsID = [trials.task_idnumber];
%    transitions = trialsID(1:end-1)-trialsID(2:end)
%    if selectID==1
%       
%    end
% end

% Initialize the truth table that will be used to further select the trials.
n_trials    = length(trials);                 % Get the remaining number of trials.
truth_table = zeros(n_trials, length(varargin)/2);
ValidTrials = zeros(n_trials, 1);

% Loop through each ..,'argument',[value],... pair on the varargin input.
for field_num = 1:2:length(varargin)
   % Check that the trial data contains the fieldname.
   if ~isempty(strmatch(lower(varargin{field_num}),lower(fields)))

    trial_elements = [trials.(varargin{field_num})];
      % Get the trial indices that have non-empty values for the desired argument.
      for n = 1:n_trials
         ValidTrials(n) = logical( ~isempty( trials(n).(varargin{field_num}) ));
      end

      % Get the desired value or range. 
      desired_range = varargin{field_num+1};
 
      if length(desired_range) == 1 % If its a scalar, make it a vector. 
         desired_range(2) =  desired_range(1);
      end
      if desired_range(1)<= desired_range(2)
         % Test each trial agains the desired value or range. 
         truth_table(logical(ValidTrials),(field_num+1)/2) = (trial_elements>=desired_range(1) & trial_elements<=desired_range(2));
      elseif desired_range(1)> desired_range(2)
         % Test each trial agains the desired value or range.
         truth_table(logical(ValidTrials),(field_num+1)/2) = (trial_elements>=desired_range(1) | trial_elements<=desired_range(2));
      elseif length(desired_range)>2 % Error checking: desired range must a two element vector [lowerlimit upperlimit].
         error(['Desired range must a two element vector ...''' varargin{field_num} ''',[lowerlimit upperlimit],...'])
      end

  % To select the spikes of a single unit, or to remove all the spikes.
   elseif strcmpi(varargin{field_num},'selectneuron')
      if isempty(prefDir) % Strip varargin from the prefDir argument and value 'cause they are not needed anymore. 
         warning('...''prefDir'',[0 180]... missing. Please provide neuron''s preferred direction.')
      end
      indices = strmatch('spike',fields);     % Find the names that begin with "spike".
      if isempty(varargin{field_num+1})       % If the value is empty remove all spikes (used when analyzing behavior only).
         trials = rmfield(trials,fields(indices(1:end)));
      elseif isempty(strmatch(varargin{field_num+1},fields,'exact')) % Error checking: Spike and channel must exist in the data.
         disp([mfilename '.m:: Valid neurons are:'])
         disp(fields(indices))
         error([varargin{field_num+1} ' not a valid neuron. '])
      else
         % Loop through each trial to get the spikes of the selected neuron and move them to a new field called "spikes".
         for k = 1:n_trials
            trials(k).spikes = trials(k).(varargin{field_num+1}); % Move the selected spikes to the "spikes" field.
            if ~isempty(prefDir) % Identify the trials in which motion direction was towards the neuron's preferred direction (1: pref, 0: nonp).
               trials(k).direction = trials(k).direction == prefDir;
            end
         end
         trials = rmfield(trials,fields(indices));   % Remove the non-selected spike fields.
      end
      truth_table(:,(field_num+1)/2) = ones(n_trials,1); % Put ones on the truth table to select all trials.
   elseif strcmpi(varargin{field_num},'align')
      alignEvent = varargin{field_num+1};
      % Just select all the trials. We'll do the alignment later.
      truth_table(:,(field_num+1)/2) = ones(n_trials,1); % Put ones on the truth table to select all trials.
   else
      disp([mfilename '.m:: Can''t select trials. Non-existent field "' varargin{field_num} '" in the trial data.'])
      truth_table(:,(field_num+1)/2) = ones(n_trials,1); % Put ones on the truth table to select all trials.
   end
end

% Remove non selected trials
remove_trials         = logical(sum(truth_table,2)~=size(truth_table,2));
trials(remove_trials) = [];

% Do time alignment ------------------------------------------------------------------
fields     = fieldnames(trials);   % Obtain valid field names.

if ~isempty(alignEvent)
   spkIndices = strmatch('spike',fields);    % Find the names that begin with "spike".
   for n = 1:length(trials)                  % Subtracts the alignTime to the timestamps of each trial.
      if ~isempty(trials(n).(alignEvent))
         alignTime = [trials(n).(alignEvent)];
         alignTime = alignTime(end);
         if isfield(trials,'eyefp_on'),  trials(n).eyefp_on       = trials(n).eyefp_on-alignTime; end
         if isfield(trials,'handfp_on'), trials(n).handfp_on      = trials(n).handfp_on-alignTime; end
         if isfield(trials,'handtouchesfix'), trials(n).handtouchesfix = trials(n).handtouchesfix-alignTime; end
         if isfield(trials,'target_dim'), trials(n).target_dim     = trials(n).target_dim-alignTime; end
         trials(n).targets_on     = trials(n).targets_on-alignTime;
         trials(n).dots_on        = trials(n).dots_on-alignTime;
         trials(n).dots_off       = trials(n).dots_off-alignTime;
         trials(n).eyefp_off      = trials(n).eyefp_off-alignTime;
         trials(n).handfp_off     = trials(n).handfp_off-alignTime;
         trials(n).saccadeini     = trials(n).saccadeini-alignTime;
         trials(n).saccadeend     = trials(n).saccadeend-alignTime;
         trials(n).reachini       = trials(n).reachini-alignTime;
         trials(n).reachend       = trials(n).reachend-alignTime;
         for i=1:length(spkIndices)          % Align the spike times too.
            trials(n).(fields{spkIndices(i)}) = trials(n).(fields{spkIndices(i)})-alignTime;
         end
      else
         RemoveTrials(n) = true;             % Will remove trials that do not have the align event.
      end
   end
end
trials(logical(RemoveTrials)) = [];             % Remove trials that do not have the align event. 
