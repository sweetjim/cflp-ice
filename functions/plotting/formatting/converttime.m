function tout = converttime(t,old_str,new_str)
t_table = [1 60 60^2 60^2*24].*[1 1/60 1/60^2 1/(60^2*24)]';
names   = {'seconds','minutes','hours','days'};
T_TABLE = array2table(t_table,'RowNames',names,'VariableNames',names);
table_idx = eval(strcat('T_TABLE.',old_str));
tout    = t*table_idx(find(strcmp(names,new_str)));
end

