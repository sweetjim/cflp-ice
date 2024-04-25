function TBOUT = Tb4thOrder(rho,S,T,varargin)
%% Interface temperature, 4th-order polynomial
%
% Finds solutions to the liquidus-temperature sensitivity balance where
%
% Sb = aTb + bTb^2 + cTb^3
%
% a = -18.7, b = -0.519, c = -0.00535;
%
% Tb is found as
%
% cTb^4+Tb^3(b+c*FqT)+Tb^2(a+b*FqT)+Tb*FqT + K = 0;
%
% where
%
% FqT   = (T.w*qw+T.i*qi+qL)./qwi;
% qi    = rho.i.*c.i*Gamma.Ti;
% qw    = rho.w.*c.w*Gamma.T;
% qL    = rho.w.*Gamma.S*property('L');
% qwi   = qi+qw;
% K     = S.w.*qL./qwi;
%%

plotting = false;
parseInput(varargin);

TBOUT = zeros(length(rho.w));

RHOW    = rho.w';
Gamma   = gammas('Ti',T.i);

c       = struct(...
    'i',property('ice','cp','T',T.i),...
    'w',property('water','cp','T',T.w));

A = -18.7;
B = -0.519;
C = -0.00535;


qi    = rho.i.*c.i*Gamma.Ti;
qw    = RHOW.*c.w*Gamma.T;
qL    = RHOW.*Gamma.S*property('L');
qwi   = qi+qw;
FqT   = (T.w*qw+T.i*qi+qL)./qwi;
K     = S.w.*qL./qwi;


Tb = linspace(T.i,0);
fn = @(Tb) C*Tb.^4+Tb.^3.*(B+C.*FqT)+Tb.^2.*(A+B.*FqT)+Tb.*FqT + K;
soln = fn(Tb);


switch size(soln,1)
    case 1
        idx = find(soln>0,1,'first');
        TBOUT = Tb(idx);
        
        if plotting
            plot(Tb,soln,Tb(idx),soln(idx),'o')
            line(Tb(idx).*[1 1],[min(ylim) soln(idx)])
            grid on
        end
    otherwise
        if plotting
            imagesc(Tb,RHOW,soln)
            addColorbar(gca,'cmap','balance','pivot',0,'latex','levels',20)
            hold on
            axis xy
            %%
            try
                contour(Tb,RHOW,soln,'k','ShowText','on')
                addlabels(gca,'latex','fs',15,...
                    'x','$T_b$ ($^\circ$C)','y','$\rho_w$ (kg m$^{-3}$)')
                set(gca,'LineWidth',1.5)
                %     latexformat
                %     print -dpng tbsoln
                shg
                axis xy
            catch TBOUT = 0;
            end
        end
        clear TBOUT
        for i=1:length(RHOW)
            TBOUT(i) =Tb(find(soln(i,:)>0,1,'first'));
        end
        hold off
end

%     
%     colorbar
%     cmocean('balance','pivot',0)
%     grid on
%%
%     its     = 10;
%     thres   = -1000;
%     map     = fn(Tb);
%     TbNZ    = Tb(map>thres&map<-thres);
%     TB      = zeros(its,2);
%     TB(1,:) = [min(TbNZ) max(TbNZ)];
%
%     for i=1:its-1
%         thres   = thres/i;
%         Tb      = linspace(min(TbNZ),max(TbNZ));
%         map     = fn(Tb);
%         TbNZ    = Tb(map>thres&map<-thres);
%
%         if ~isempty(TbNZ)
%             TB(i,:) = [min(TbNZ) max(TbNZ)];
%         else
%             TB(i,:) = NaN.*[1 1];
%         end
%         if isempty(min(TbNZ))||isempty(max(TbNZ))
%             TB(i:end,:) = [];
%             break
%         end
%     end
%    TBOUT(j) = mean(TB(end,:));
%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'plot','plotting'}
                    plotting = true;
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
