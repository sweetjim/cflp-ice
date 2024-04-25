function [b,varargout] = iw_iterate(varargin)
%% IW_ITERATE 
%   H:      domain depth / height
%   L:      domain bounds [min max]
%   res:    resolution (m)
%   dt:     time resolution (s)
%   h:      topography height (m)
%   k0:     fundamental wavenumber of h (1/m)
%   N:      stratification (rad/s)
%   u:      linear velocity (m/s)
%   ut:     oscillating velocity magnitude (m/s)
%   omega:  oscillation frequency (rad/s)
%   f:      Coriolis frequency (rad/s)
%   nu:     viscosity (m2/s)
%   n:      number of harmonics to compute
%   alpha0: harmonic regime
%   hfcn:   custom wavemaker profile
%% Input
global H res omega N u ut nu n k0 alpha0 f dt xres zres
%#ok<*GVMIS> 
iterate = true;
outputs = {};
H       = 1;
L       = [-1 1]*5;
k0      = sqrt(diff(L)/(max(L)^.5));
h0      = H/10;
res     = 1e-3;
alpha0  = 2;
f       = 0; 
dt      = .3;
N       = .7;
nu      = 1e-5;
n       = 2;
omega   = N/2;
u       = 0;
ut      = 20e-3;
x0      = 0;
xres    = [];
zres    = [];
hfcn    = [];
parseInput(varargin)
n       = 0:max(n);
z       = linspace(0,H,1/res);
x       = linspace(L(1),L(2),1/res);
if ~isempty(xres)
x       = linspace(L(1),L(2),xres);
end
if ~isempty(zres)
z       = linspace(0,H,zres);
end

h_fnc   = @(h,efolding,x0) h*exp(-((x-x0)).^2/(efolding)^2);
h       = h_fnc(h0,1/k0,x0);
if ~isempty(hfcn)
h       = hfcn;
end
T       = (2*pi/omega);

%% Printing
alpha = omega/sqrt(N^2-omega^2);
theta = atan(alpha);
% fprintf('Theta = %.2f degrees\n',rad2deg(theta))
%% Solver
iw          = iw_initialize(H,z,x,h,N,u,ut,omega,f,nu,n,alpha0);
iw.N        = N;
iw.x        = x;
iw.omega    = omega;
bn = iw.bnn;
b = iw.b';
w = iw.w';
u = iw.u';

if nargout==0
    clf
    imagesc(x,z,iw.b')
    hold on
    fill(x,H-h,'k')
    hold off
    axis ij
end

tvec            = linspace(0,T,ceil(T/dt));
B = zeros(size(b,1),size(b,2),2);
B(:,:,1) = b;
B(:,:,2) = timestepper(iw,2);

% dwdt = diff(B,1,3)/mean(diff(tvec));
% % [dwdx,dwdz,dwdt]=gradient(B,x,z,tvec(1:2));
% eta         = dwdt;%#ok<*NASGU> %b/N^2;%; % displacement field
% amplitude   = max(-u/(N*sin(atan(omega/N)))/omega,[],'all');

if isempty(outputs)
    varargout{1} = w;
    return
end

if iterate
    iw.tvec = tvec;
    b_time_field    = zeros(numel(z),numel(x),numel(tvec));

    if ~isempty(gcp('nocreate'))
        b_time_field = permute(b_time_field,[3 1 2]);
        q = parallel.pool.DataQueue;
        afterEach(q,@displayProgressparfor)
        displayProgressparfor([],numel(tvec))
        loop = length(tvec);
        iw.tvec = tvec;
        parfor t=1:loop
            b_time_field(t,:,:) = iw_timestep(iw,t);%(:,frame(1):frame(2));      % save the 'lab' frame
            send(q,[])
        end
        b_time_field = permute(b_time_field,[2 3 1]);
    else
        for t=1:length(tvec)
            b_time_field(:,:,t) = iw_timestep(iw,t);%(:,frame(1):frame(2));      % save the 'lab' frame
            displayProgress('Iterating',t,1,numel(tvec))
        end
    end
    b = b_time_field;
    t = tvec;
    [~,~,eta] = gradient(b,x,z,t);
end



for ii=1:numel(outputs)
    exec = sprintf('varargout{%i}=%s;',ii,outputs{ii});
    eval(exec)
end
%% Functions
    function b_t = timestepper(iw,t)
        i = complex(0,1);
        %% Buoyancy
        b0      = iw.w0.*N^2./(i.*(iw.OMEGA0));
        b1      = iw.w1.*N^2./(i.*(iw.OMEGA1));
        bhat    = (b0.*iw.what0.*exp(i*iw.OMEGA0*tvec(t))+b1.*iw.what1.*exp(i*iw.OMEGA1*tvec(t)));
        
        bhat(iw.nn==0)=bhat(iw.nn==0)/2;                              % remove duplicate n=zeros
        bhat(isnan(bhat))       = 0;
        bhat(isinf(abs(bhat)))  = 0;
        
        bb              = iffts(bhat,length(x),1,length(x));    % back to real space
        b_t             = squeeze(sum(bb(:,:,:),3))';           % sum over harmonics
    end
%% Input parser
    
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'H','res','omega','N','u',...
                        'ut','nu','n','k0','h0','alpha0','f','dt',...
                        'x0','xres','zres'}
                    exec = sprintf('%s=namevalue;',items{m});
                    evalin('caller',exec)
                case 'hfcn'
                    hfcn = namevalue;
%                 case 'xres'
%                     xres = namevalue;
%                 case 'zres'
%                     zres = namevalue;
                case 'L'
                    L = namevalue;
                    if numel(L)==1
                       L = [0 L]; 
                    end
                case {'out','output'}
                    outputs = items(m+1:end);
                    return
                case {'initial','t0','zero'}
                    iterate = false;
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

