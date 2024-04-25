function Root = FixedPointIteration(fh,x1,varargin)
%   This function uses fixed point iteration to find the root of a 
%   nonlinear function given the following set of criteria. 
%   - The function magnitude remains always less than or equal to 1
%   - The function can be written in the form x = f(x), such as the
%     case of fh = @(x) x-cos(x), one could write x = cos(x). The fixed
%     point iteration would then ouput the x such that this condition is
%     satisfied, which will be the root.
%
% Required input arguments:
%
%    fh:        Function handle, either a scripted function or an anonymous
%               function. An example of an anonymous function is given here
%
%                           fh = @(x) 5*sqrt(x) - x
%
%    x1:        The first guess for the solver
%               
%  Optional input arguments:
%
%  iteration:   The maximum number of attempts that the solver will 
%               make in searching. This is fast, but can easily diverge, so
%               limited iterations are recommended.
%  tolerance:   The accuracy that the solution must have.
%
if isempty(varargin)
    iterations = 100;
    tol = 1e5*eps;
elseif length(varargin{:})==1
    iterations = varargin{1};
    tol = 1e5*eps;
else
    iterations = varargin{1};
    tol = varargin{2};
end
x2 = fh(x1);
for i = 1:iterations
    if abs(x2-x1)>tol
        x1 = x2;
        x2 = fh(x1);
    else
        Root = x2;
        break
    end
end
end