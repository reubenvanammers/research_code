function [Time, strain, restoring] =  strain_restoring(fhandle,varargin)
global external_force restoring_rec restoring_t_rec
[Time,Y]=fhandle(varargin{:});
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = max(xvalues,[],2)/(max(xvalues(1,:))-min(xvalues(1,:)));
restoring_temp = restoring_rec./external_force;
[restoring_t_rec,ia,ic] = unique(restoring_t_rec);
restoring_temp = restoring_temp(ia);%deletes duplicate time entries for interpolation
restoring = interp1(restoring_t_rec,restoring_temp,Time);%interpolates restoring force to be same size as time vector