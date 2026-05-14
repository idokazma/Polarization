function out = merge_defaults(defaults, overrides)
% merge_defaults  Shallow-overlay overrides onto defaults.
% Every field in overrides replaces the same-named field in defaults; all
% other fields of defaults pass through untouched.
%
% Usage:
%   cfg = merge_defaults(defaults, overrides);
%
% If overrides is empty or missing, defaults is returned unchanged.

out = defaults;
if nargin < 2 || isempty(overrides), return; end

fn = fieldnames(overrides);
for k = 1:numel(fn)
    out.(fn{k}) = overrides.(fn{k});
end

end
