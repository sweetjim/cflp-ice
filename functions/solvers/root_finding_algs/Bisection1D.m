function Root = Bisection1D(fh,Low,High,varargin)
%   This function uses the Bisection Method to find the root of a 
%   nonlinear function given the following set of criteria. 
%   - The bounds are opposite signs, or one is equal to zero. This 
%     function will error when given bounds that are  equal in sign, as 
%     this could imply multiple roots or none at all
%   - The bounds will preferably not contain discontinuities, as the 
%     direction that the bisection method should take will be unclear
%
% Required input arguments:
%
%    fh:        Function handle, either a scripted function or an anonymous
%               function. An example of an anonymous function is given here
%
%                           fh = @(x) 5*sqrt(x) - x
%
%   Low:        The lower bound in which to look for the root
%
%  High:        The upper bound in which to look for the root
%               
%  Optional input arguments:
%
%  tolerance:   The accuracy that the solution must have, shown as the
%               difference of the calculated midpoint from zero
%
a = fh(High);
b = fh(Low);
if isempty(varargin)
    tol = 1e5*eps;
else
    tol = varargin{1};
end
if sign(b)==sign(a)
    disp(['These bounds contain either no roots, a discontinuity,' ...
    ' or multiple roots'])
    disp('Please adjust bounds')
    Root = nan;
    return
end
Mid = (Low+High)/2;
c = fh(Mid);
if isinf(c)||isnan(c)
    disp('Discontinuity in the within the boundaries')
    Root = nan;
    return
end
while abs(c)>tol
    if sign(c)~=sign(a)
        Low = Mid;
        Mid = (Mid+High)/2;
    else
        High = Mid;
        Mid = (Mid+Low)/2;
    end
    c = fh(Mid);
    if isinf(c)||isnan(c)
        disp('Failed to converge')
        Root = nan;
        return
    end
end
Root = Mid;
end