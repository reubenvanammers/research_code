function [Time, strain, restoring] =  strain_restoring(fhandle,varargin)
%calculates the strain and restoring force for a given model fhandle with
%parameters varargin. Model needs to have following global variables baked
%in in order to work. 
global restoring_rec t_rec

[Time,Y]=fhandle(varargin{:});
l = length(t_rec)-length(restoring_rec);
t_rec = t_rec(l+1:end);
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
restoring_temp = restoring_rec;
[t_rec,ia,~] = unique(t_rec);
restoring_temp = restoring_temp(ia);%deletes duplicate time entries for interpolation
restoring = interp1(t_rec,restoring_temp,Time);%interpolates restoring force to be same size as time vector