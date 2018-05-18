function [value,isterminal,direction] = strain_event(t,y)
%Event checks whether to stop stop applying strain, as system has
%equilibriated

global stress_rec t_rec

if t_rec(end)-t_rec(end-1) ~= 0
    stress_gradient = (stress_rec(end)-stress_rec(end-1))/(t_rec(end)-t_rec(end-1));
    has_equilibriated = abs(stress_gradient)<1e-10;
else
    has_equilibriated = false;
end

value = double(~has_equilibriated);
isterminal=1;
direction = 0;
