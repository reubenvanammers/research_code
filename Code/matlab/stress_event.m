function [value,isterminal,direction] = stress_event(t,y)
%Event checks whether to stop stop applying stress, as system has
%equilibriated
global external_force restoring_rec strainflag maxlength
%If restoring force sufficiently close to the applied force, equilibriated,
%so stop the system
has_equilibriated = (abs(restoring_rec(end)-external_force))<(1e-10/external_force);

%If there is a flag to see whether cells have gone above certain strain,
%call for system to stop


if ~isempty(maxlength)
    [V,~] = matricize(y);
    Vx = V(:,1);
    if max(Vx)-min(Vx) > maxlength
        strainflag = true;
    end
end

if ~isempty(strainflag)
    has_broken = strainflag;
else
    has_broken = false;
end

value = double(~(has_equilibriated || has_broken));
isterminal = 1;
direction = 0;