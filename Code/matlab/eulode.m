function [t,y] = eulode(dydt, t, y0, asdf)



n = length(t);
%if necessary, add an additional value of t 
%so that range goes from t=ti to tf

y = ones(n,length(y0));
y(1,:) = y0;

for i = 1:n-1 %implement Euler's Method

    y(i+1,:) = y(i,:) + dydt(t(i),y(i,:))'*(t(i+1)-t(i));
end