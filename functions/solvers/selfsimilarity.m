function alpha = selfsimilarity(y1,t1,y2,t2)
%SELFSIMILARITY 
% A self similar process is one that satisfies:
%   y(t) = a^alpha*y(t/a),
% where "a" is a rescaling factor and "alpha" is the self-similarity
% parameter.
% The self-similarity parameter may be determined by the test:
%   alpha = [ln(s2)-ln(s1)]/[ln(n1)-ln(n2)]
% where "s" denotes the standard deviation of the probability distribution
% and "n" denotes the size of the independent variable
% 
%% Window algorithm
% ???



%%
tiledlayout(1,2)
n(1)=nexttile;
plot(n(1),t1,y1,t2,y2)
n(2)=nexttile;
bins = 50;
s1 = std(histcounts(y1,bins,'Normalization','probability'));
s2 = std(histcounts(y2,bins,'Normalization','probability').Values);

n1 = diff([min(t1) max(t1)]);
n2 = diff([min(t2) max(t2)]);

alpha = (log(n2)-log(n1))/(log(s2)-log(s1));
plot(n(1),t1,y1,t2,y2,t1,y2*1^alpha)

%%
loglog([n1 n2],[s1 s2])
    function windowSampling(t1,t2)

    end
end

