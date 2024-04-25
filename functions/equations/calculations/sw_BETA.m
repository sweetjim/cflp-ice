function out = sw_BETA(S,T)
%% Seawater haline contraction coefficient beta
% Calculates beta at atmospheric pressure from the seawater package but
% for non-equal S,T vector sizes.
%%
if numel(S)>1&&numel(T)==1
    %% S-vector
    for i=1:length(S)
       out(i) = sw_alpha(S(i),T,0); 
    end
elseif numel(S)==1&&numel(T)>1
    %% T-vector
    for i=1:length(T)
       out(i) = sw_alpha(S,T(i),0); 
    end
elseif numel(S)==1&&numel(T)==1
    %% Value
    out     = sw_beta(S,T,0);
else
    %% ST-matrix
    if size(S,1)==size(S,2)&&size(T,1)==size(T,2)
        s = S;
        t = T;
    else
        [s,t]   = meshgrid(S,T);
    end
    out     = sw_beta(s,t,0);
end

