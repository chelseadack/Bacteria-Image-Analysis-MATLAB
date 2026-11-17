function positive_threshold = func_gate_threshold(negative_values, gate_percent)
% FUNC_GATE_THRESHOLD  GFP positive gate set from the negative control cells.
%
%   positive_threshold = func_gate_threshold(negative_values, gate_percent)
%
%   A cell is GFP positive if it is brighter than all but gate_percent % of the negative control cells,
%   eg. gate_percent = 5 sets the gate at the 95th percentile of the negative control.
%
%   Input parameters:
%       - negative_values = fluorescence of every negative control cell.
%       - gate_percent = % of negative control cells allowed above the gate (5 or 1 are typical).
%   Output:
%       - positive_threshold = cells with fluorescence above this value are GFP positive.

    arguments
        negative_values {mustBeNumeric}
        gate_percent (1,1) {mustBeNumeric, mustBeInRange(gate_percent, 0, 100, 'exclude-lower')}
    end
    if isempty(negative_values)
        error('func_gate_threshold:noCells', 'No negative control cells to set the gate from.');
    end

    negControl = sort(negative_values(:));
    positive_threshold = negControl(ceil((1 - gate_percent / 100) * end));

end
