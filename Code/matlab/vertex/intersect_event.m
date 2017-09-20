function [value,isterminal,direction] = intersect_event(t,y)

global C
[V,V_ref] = matricize(y);
value = is_self_intersecting(C,V) || is_self_intersecting(C,V_ref);

value = double(value);
direction = 0;
isterminal = 1;