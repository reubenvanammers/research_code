function included_cell = cell_inclusion(V,C)
%Creates lists of vertices that are included in each cell. Used in internal vertex calculation, and is precalculated
%for speed
included_cell = cell(1,length(V));
m = max(cellfun(@length,C)); 
for l = 1:length(C);
    while length(C{l})<m
        C{l} = [C{l} NaN];
    end 
end
for i = 1:length(V);
    [cells,~] = find(cell2mat(C)==i);
    included_cell{i} = cells';
end