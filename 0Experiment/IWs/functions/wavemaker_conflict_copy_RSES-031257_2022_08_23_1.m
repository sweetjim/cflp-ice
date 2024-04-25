function wavemaker(N,theta,amplitude_mm)
%% wavemaker generates waveforms for use with the wavemaker program
% 
%% Assumed geometry
k0 = 2*pi/sqrt(5e-2);

%% Input
% theta = @(N,omega) rad2deg(asin(omega/N));
omega   = N*sin(deg2rad(theta));
t_res   = 50e-3; % 50ms precision
period  = 2*pi/omega;
t       = linspace(0,round(period,3),period/t_res);
y       = sin(omega*t)*amplitude_mm/(24.5/2);
y(end)  = 0;
plot(t,y)
shg
%% 
fr = (amplitude_mm*1e-3)*k0/N;
fprintf('Froude number: %.2f\n',fr)
return
%% Output
writetable(...
    table(y'),...
    strcat(strrep(sprintf('N_%.2frads-theta_%.2fdeg-amp_%.1fmm-Fr_%.2f',...
    N,theta,amplitude_mm),'.','_'),'.csv'),...
    'WriteVariableNames',false)
end

