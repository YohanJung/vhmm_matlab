function [row,col]=find_arcs(True_inter)

p = size(True_inter,1);
row = []; col = [];

for i = 1:p
    for j = 1:p
        if  True_inter(i,j)
            row = [row;i];
            col = [col;j];
        end
    end
end
