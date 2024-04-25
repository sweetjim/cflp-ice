function OUT = staircase(z,rho,varargin)
%% Staircase
% Converts a linear density field into a staircase density field.
%
% Parameters: z [vec], rho [vec], steps [int or double (for h)]
% 
%  example:
%   z = linspace(0,1);rho = linspace(1000,1050); 
%   staircase(z,rho,'h','fixh','step',-.08:.02:1);
%%

fixh    = true;
type    = 'n';
n       = 10;
parseInput(varargin);

if fixh
   n0 = cumsum(abs(n));
   n  = n(n0<1);
end

switch type
    case 'n'
        try
            X   = rho;
            rho = [rho(1) rho(end)];
            rho = linspace(rho(1),rho(2),n);
            OUT = interp1(rho(1):diff(rho):rho(end),rho,X,'next');
        catch
            error('n is expected to be an integer.')
        end
    case 'h'
        if length(n)>1
            %% h(z)
            n = abs(n);
            if sum(n)>(z(end)-z(1))
                error('total h cannot exceed the bounds of z')
            end
            
            ni = cumsum(n);
            clear zidx
            zidx(1) = 1;
            for i=1:length(ni)
                if isempty(find(z>ni(i),1,'first'))
                    zidx(i+1) = length(rho);
                else
                    zidx(i+1) = (find(z>ni(i),1,'first'));
                end
                
                if i>1
                    rho(zidx(i):zidx(i+1)) = rho(zidx(i+1));
                end
            end
            rho(1:zidx(2)) = rho(1);
            rho(zidx(end):end) = rho(end);
            
            OUT = rho;
            
            
        else
            %% constant h
            n = abs(n);
            n0 = (z(end)-z(1))/n;
            
            if n0==2
               staircase(z,rho,3)
               return
            end
            
            if mod(1,n0)==0 || mod(n0,1)==1
                X   = rho;
                rho = [rho(1) rho(end)];
                rho = linspace(rho(1),rho(2),n0);
                OUT = interp1(rho(1):diff(rho):rho(end),rho,X,'next');
            else
                n0 = 0;
                i = 1;
                while cumsum(n0)<=1
                    n0(i+1) = n0(i)+n;
                    i = i+1;
                end
                
                zidx(1) = 1;
                for i=1:length(n0)
                    zidx(i+1) = (find(z>n0(i),1,'first'));
                    if i>1
                        rho(zidx(i):zidx(i+1)) = rho(zidx(i+1));
                    end
                end
                rho(1:zidx(2)) = rho(1);
                rho(zidx(end):end) = rho(end);
                
                OUT = rho;
            end
        end
end


plot(z,OUT);
grid on

%%
    function parseInput(varargin)
        %%
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'n'}
                    type = 'step';
                case {'h'}
                    type = 'h';
                case {'fixh'}
                    fixh = true;   
                case {'step'}
                    n = namevalue;
            end
            m = m+1;
            if m>length(items);break;end
        end
        %%
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

