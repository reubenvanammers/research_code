function stat =connected_cells(C,V)
N = size(V,1);
G = zeros(N);
stat = true;
for i = 1:length(C);
    for k = 1:length(C{i})
        v1 = C{i}(k);
        v2 = C{i}(mod(k,length(C{i}))+1);
        G(v1,v2) = 1;
        G(v2,v1) =1;
    end
end
G =sparse(G);
[S,C]=graphconncomp(G);

stat = S==1;