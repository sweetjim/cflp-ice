% WSSIZE
%
% Display size of variables in the base workspace using B, KB, MB and GB
function wssize(variable)
if ~nargin
    % Obtain base workspace contents
    ws_contents = evalin('base', 'whos');
    usevariable = false;
else
    if isempty(variable)
        warning('Empty variable')
        return
    end
    ws_contents = whos('variable');
    usevariable = true;
    validclasses = {'struct','icedata','iwdata'};
end
% Loop through contents and display size on screen
[~,sorting]=sort([ws_contents.bytes]);
ws_contents = ws_contents(sorting);
for i = 1:length(ws_contents)
    if usevariable
        if contains(ws_contents.class,validclasses)
            fields = fieldnames(variable);
            dims   = max(ws_contents.size);
            cur_size = zeros(size(fields,1),dims);
            for k=1:dims
                command = cellfun(@(x) sprintf('variable(%i).%s;',k,x),fields,'UniformOutput',false);
                for j=1:numel(command)
                    tmp      = evalin('caller',command{j}); %#ok<NASGU>
                    cur_size(j,k) = whos('tmp').bytes;
                end
            end
            tbl = table(fields,cur_size,'VariableNames',{'fields','size'});
            tbl = sortrows(tbl,"size",'ascend');

            %%
            b = barh(categorical(tbl.fields,tbl.fields),tbl.size/1024^2,'stacked');
            %             set(gca,'Xscale','log')
            xtips1 = b(end).YEndPoints;
            ytips1 = b(end).XEndPoints;

            labels1 = arrayfun(@(x) sprintf('%.1f MB',x),sum(tbl.size,2)/1024^2,'UniformOutput',false);
            text(xtips1,ytips1,labels1,'VerticalAlignment','middle')
            xlabel('Size (MB)')
        end
    end
    cur_size = ws_contents(i).bytes;
    if cur_size > 1024^3
        fprintf('%-15s: %8.3f GB\n', ws_contents(i).name, cur_size/1024^3);
    elseif cur_size > 1024^2
        fprintf('%-15s: %8.3f MB\n', ws_contents(i).name, cur_size/1024^2);
    elseif cur_size > 1024
        fprintf('%-15s: %8.3f KB\n', ws_contents(i).name, cur_size/1024);
    else
        fprintf('%-15s: %8.3f  B\n', ws_contents(i).name, cur_size);
    end
end
end