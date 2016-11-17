function circularity_vec = circularity(C,V)

M = length(C);
real_cell_areas = zeros(1,M);
real_cell_circumferences = zeros(1,M);

for i = 1:M
    real_cell_areas(i) = cell_area(i,C,V);
    real_cell_circumferences(i) = cell_circumference(i,C,V);
end

circularity_vec = 4*pi*real_cell_areas.^(2).*real_cell_circumferences.^(-2);