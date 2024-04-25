function out = density_find(var_str,alt_var_tgt,density_tgt)
%% density_find calculates the value of the variable (S OR T)
% required to achieve the target density.
% Parameters:
%   var_str     'T' or 'S'
%   alt_var_tgt  Alternate variable target value    
%   denisty_tgt  Target denisty
%
% Example:
%  Find the salinity at 10 degC and 1010 kg/m3
%   density_find('S',10,1010)
%   ans = 13.2171 (psu)
%  Check
%   density(density_find('S',10,1010),10)
%   ans = 1010 (kg/m3)
%%

switch var_str
    case 'T'
        field = linspace(-5,50,1e5);
        f = density(alt_var_tgt,field);
        g = f>density_tgt;
        if numel(find(g))>1
            out = [field(find(g,1,'first')) field(find(g,1,'last'))];
            return
        end
        out = field(find(g,1,'first'));
    case 'S'
        field = linspace(0,200,1e5);
        f = density(field,alt_var_tgt);
        if numel(density_tgt)>1
            g = arrayfun(@(x) find(f>x,1,'first'),density_tgt,'UniformOutput',false);
            out = field(cell2mat(g));
            return
        end
        g = f>density_tgt;
        out = field(find(g,1,'first'));
end
end

