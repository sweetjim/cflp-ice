function out = density_find(var_str,alt_var_tgt,density_tgt)
%% density_find calculates the value of the variable (S OR T)
% required to achieve the target density.
%%
switch var_str
    case 'T'
        f = density(alt_var_tgt,linspace(-20,50,1e5));
        g = f>density_tgt;
        if numel(find(g))>1
            out = [f(find(g,1,'first')) f(find(g,1,'last'))]
            return
        end
    case 'S'
end
end

