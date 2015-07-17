function ppe = searchpp(pp, searchfor)
% ppe = searchpp(pp)
% search output of readpp by name

idx = false(1, length(pp));

if ~iscell(searchfor)
    searchfor = {searchfor};
end

% if numel(pp) == 1
%     fields = fieldnames(pp);
% end

for s = searchfor
    idx = idx | cellfun(@(x) ~isempty(x), regexp({pp.name}, s));
end

ppe = pp(idx);
