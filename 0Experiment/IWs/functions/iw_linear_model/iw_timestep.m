function [b_t,u_t,w_t] = iw_timestep(iw,t)
i       = complex(0,1);
N       = iw.N;
tvec    = iw.tvec;
omega   = iw.omega;
x       = iw.x;
%% Buoyancy
b0      = iw.w0.*N^2./(i.*(iw.OMEGA0));
b1      = iw.w1.*N^2./(i.*(iw.OMEGA1));
bhat    = (b0.*iw.what0.*exp(i*iw.OMEGA0*tvec(t))+b1.*iw.what1.*exp(i*iw.OMEGA1*tvec(t)));

bhat(iw.nn==0)=bhat(iw.nn==0)/2;                              % remove duplicate n=zeros
bhat(isnan(bhat))       = 0;
bhat(isinf(abs(bhat)))  = 0;

u_t   = iw.u*cos(omega*t);
w_t   = iw.w*cos(omega*t);
bb    = iffts(bhat,length(x),1,length(x));    % back to real space
b_t   = squeeze(sum(bb(:,:,:),3))';           % sum over harmonics
end