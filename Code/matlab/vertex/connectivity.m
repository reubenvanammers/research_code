function connectivitylist = connectivity(C,V)
for v = 1:length(V)
    for m = 1:length(C)
        if any(C{m}==v)
            connectivitylist(m,v) = 1;
        end
    end
end
end