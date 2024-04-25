function [Tb,Sb] = TSinterface(s,t,varargin)
%% Temperature-salinity interface conditions
% Approximates the interface properties.
%
% Parameters: Sinf, Tinf, (Optional) Ti (-5), 
% (Optional) method ('intermediate','liquidus','heatflux','both')
% (Optional) sb [scalar guess] (def. 0.5)

%%
Ti      = -5;
method  = 'liquidus';
sb      = 0.5;
parseInput(varargin);

Sb = s*sb;  % initial guess

switch method
    case 'liquidus'
        Tb = liquidus(Sb,0);
    case 'heatflux'
        Tb = Tb_heatflux(Sb,s,Ti,t);
    case 'both'
        %% Iterative process
        for i=1:20
            Tb = [liquidus(Sb,0) Tb_heatflux(Sb,s,Ti,t)];   % liquidus and heat-flux approx.
            Sb = mean(Tb);                                  % take midpoint
            TB(i,:) = Tb;
            SB(i) = Sb;
        end
        Sb = rms(SB);                                      %
        Tb = mean([liquidus(Sb,0) Tb_heatflux(Sb,s,Ti,t)]);
        
    case 'inter'
        Sb = s/2;
        Tb = t/2;
end

%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'Ti','T_i'}
                    Ti      = namevalue;
                case {'Liquidus','liquidus'}
                    method = 'liquidus';
                case {'inter','intermediate'}
                    method = 'inter';
                case {'heatflux','hf'}
                    method = 'heatflux';
                case {'both','lhf','hfl'}
                    method = 'both';
                case {'sb'}
                    sb = namevalue;
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

