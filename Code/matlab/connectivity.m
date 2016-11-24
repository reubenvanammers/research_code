function connectivitylist = connectivity(C,V)
for v = 1:length(V)
    for c = 1:length(C)
        if any(C{c}==v)
            connectivitylist(c,v) = 1;
        end
    end
end
end