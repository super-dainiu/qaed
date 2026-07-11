function tf = qaed_accel(setting)
% QAED_ACCEL Toggle dispatch to the C++ core (cpp/qaed_mex).
%
% tf = QAED_ACCEL()        returns true if the C++ core will be used.
% QAED_ACCEL(false)        force the pure-MATLAB reference implementations.
% QAED_ACCEL(true)         re-enable the C++ core (default).
% QAED_ACCEL('off'/'on')   same as above.
%
% The C++ core is used only when it is both enabled AND qaed_mex is on the
% path (built via scripts/build_mex_and_check.sh), so code always works.
persistent state
if isempty(state), state = true; end
if nargin > 0
    if ischar(setting) || isstring(setting)
        setting = ~strcmpi(char(setting), 'off');
    end
    state = logical(setting);
end
tf = state && exist('qaed_mex', 'file') == 3;
end
