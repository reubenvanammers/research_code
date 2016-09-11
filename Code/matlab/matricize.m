function [R,P] =  matricize(column)
%converts columnn of x,y data to xy coordinates in 2d matrix. Inverse of
%columnize.
N = length(column);
N = N/2;
real_column = column(1:N);
reference_column = column(N+1:2*N);
R = reshape(real_column,[],2);
P = reshape(reference_column,[],2);