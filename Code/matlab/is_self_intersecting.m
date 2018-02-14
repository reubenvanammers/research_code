%Checks to see whether or not a certain cell intersects with itself so that
%the simulation can stop.
function b = is_self_intersecting(C,V)
b = false;
for i = 1:size(C,1)
    x = V(C{i},1);
    y = V(C{i},2);
    [~,~,segments] = selfintersect(x,y);
    if ~(isequal(segments,double.empty(0,2))||isequal(segments,[]))
        b = true;
        segments
        i
    end
end