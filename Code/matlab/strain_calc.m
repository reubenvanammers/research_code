function [Time,strain,Y] = strain_calc(fnhandle,varargin)
%calculatse strain for fnhandle with arguments varargin.
[Time,Y]=fnhandle(varargin{:});
N = size(Y,2)/4;
xvalues = Y(:,1:N);
strain = (max(xvalues,[],2)-min(xvalues(1,:)))/(max(xvalues(1,:))-min(xvalues(1,:)));
figure
plot(Time,strain,'r')
