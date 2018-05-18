function [value,isterminal,direction] = intersect_event(t,y)
%creates event if any cells self intersect. Usually stops code if this is
%the case if used as an event in the code
global C
[V,V_ref] = matricize(y);
value = is_self_intersecting(C,V) || is_self_intersecting(C,V_ref);

value = double(value);
direction = 0;
isterminal = 1;