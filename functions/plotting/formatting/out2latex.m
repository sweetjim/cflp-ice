function out2latex(folderpath,filename,varargin)
%% Out to latex
% Formats figures to be suitable for latex documents.
% 
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   filename    :   [char]  - name of figure output 
% 
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%   'path'      :   [char]  - folder path of output (default is cwd)
%   'type'      :   [char]  - type of file output
%       'png' (optional transparent background), 'jpeg', 'tiff'
%   'transparent'   [char]  - enable transparency 
%   'resolution':   [double]- set the resolution in ppi 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
% type    = '-dpng';
type    = 'png';
fig     = gcf;
res     = 864;
istransparent = false;
isfulltransparent = false;


parseInput(varargin)

% if ~strcmp(type,'-dpng')
%    type = strcat('-d',type); 
% end
% if ~isempty(folderpath)
filename =  fullfile(folderpath,filename);
% end

if istransparent||isfulltransparent
    set(fig,'Color','none')
    if isfulltransparent
       set(findobj(fig,'Type','Axes'),'Color','none') 
    end
    export_fig(fig,filename,'-png','-transparent',sprintf('-r%i',floor(res)))
    set(fig,'Color','w')
    set(findobj(fig,'Type','Axes'),'Color','none') 
    return
end
warning off
export_fig(fig,filename,strcat('-',type),sprintf('-r%i',floor(res)))
warning on
% pause(.5)
% print(fig,filename,type)


%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
%                 case {'path','folder','root'}
%                     folderpath = namevalue;
                case {'type','file'}
                    type = namevalue;
                case 'transparent'
                    istransparent = true;
                case 'fulltransparent'
                    isfulltransparent = true;
                case isa(items{m},class(gcf))
                    fig = items{m};
                case {'r','resolution','res'}
                    res = namevalue;
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

