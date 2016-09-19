function [C,V,V2] = t1swap(vertex1,vertex2,dsep,C,V,varargin)
%Executes t1swap between vertex1 and vertex2. Can swap between reference
%cells as well if added as an extra argument.


%Calculates new position of cells after t1swap.
avpos = 0.5*(V(vertex1,:)+V(vertex2,:));
seperation_vector = V(vertex1,:)-V(vertex2,:);
new_seperation = [0 1;-1,0]*seperation_vector'./norm(seperation_vector);
V(vertex1,:) = avpos+dsep.*new_seperation';
V(vertex2,:) = avpos-dsep.*new_seperation';

if nargout ==3 %ie if want to swap reference cells as well
    V2 = varargin{:};
    avpos = 0.5*(V2(vertex1,:)+V2(vertex2,:));
    seperation_vector = V2(vertex1,:)-V2(vertex2,:);
    new_seperation = [0 1;-1,0]*seperation_vector'./norm(seperation_vector);
    V2(vertex1,:) = avpos+dsep.*new_seperation';
    V2(vertex2,:) = avpos-dsep.*new_seperation';
end

[including_cells_1,~] = find(cell2mat(C)==vertex1);
[including_cells_2,~] = find(cell2mat(C)==vertex2);
common_cells = intersect(including_cells_1,including_cells_2);
adjacent_cells = setdiff(union(including_cells_1,including_cells_2),common_cells);
no_of_cells = length(common_cells)+length(adjacent_cells);
total_cells = reshape([common_cells';adjacent_cells'],1,[]);%cells A,B,C,D in Fletch2013Implementing respectively


%Restructures C so that the new connectivity is applied.
for i = 1:no_of_cells;
    if i==1 || i==3; %A or C
        l = total_cells(i);
        len = length(C{l});
        k1=find(C{l}==vertex1)-1;
        n11 = C{l}(mod(k1 - 1,len)+1);
        n12 = C{l}(mod(k1 + 1,len)+1);
        if norm([norm(V(vertex1,:)-V(n11,:)),norm(V(vertex1,:)-V(n12,:))])< norm([norm(V(vertex2,:)-V(n11,:)),norm(V(vertex2,:)-V(n12,:))])
            C{l} = setdiff(C{l},vertex2,'stable');
        else
            C{l} = setdiff(C{l},vertex1,'stable');
        end
    else
        %B or D
        l = total_cells(i);
        len = length(C{l});
        internal_vertex = intersect([vertex1 vertex2],C{l});
        external_vertex = setdiff([vertex1 vertex2],internal_vertex);
        k1=find(C{l}==internal_vertex)-1;
        n11 = C{l}(mod(k1 - 1,len)+1);
        if norm(V(n11,:)-V(internal_vertex,:)) > norm(V(n11,:)-V(external_vertex,:))
            if k1==0;
                tempcell = [external_vertex C{l}];
            else
                tempcell = [C{l}(1:k1) external_vertex C{l}((k1+1):end)];
            end
        else
            if k1==len-1
                tempcell = [C{l} external_vertex];
            else
                tempcell = [C{l}(1:(k1+1)) external_vertex C{l}((k1+2):end)];
            end
        end
        C{l} =tempcell;
    end
end

               