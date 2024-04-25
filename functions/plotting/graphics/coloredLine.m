function s = coloredLine(X,Y,C,ax)
%COLOREDLINE plot a coloured lineplot of vector Y versus vector X with colour C.
% Parameters: (Required)
%   X   -   X vector [1xn column vector]
%   Y   -   Y vector [1xn column vector]
%   C   -   C vector [1xn column vector]
% Parameters: (Optional)
%   ax  -   Parent axes (default is gca)
%
% Outputs: (Optional)
%   s   -   surface data
%%
if nargin<4
    ax = gca;
end
X = getColumn(X);
Y = getColumn(Y);
C = getColumn(C);
x = zeros(size(X));
s=surface(ax,[X X],[Y Y],[x x],[C C], ...
    'FaceColor','none', ...
    'EdgeColor','interp', ...
    'LineWidth',2);
if ~nargout,clear s,end
    function in = getColumn(in)
        if isrow(in)
            in = in';
        end

    end
end

