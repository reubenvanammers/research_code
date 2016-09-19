function included_cell = cell_inclusion(V,C)
included_cell = cell(1,length(V));
for i = 1:length(V);
    [cells,~] = find(cell2mat(C)==i);
    included_cell{i} = cells';
end