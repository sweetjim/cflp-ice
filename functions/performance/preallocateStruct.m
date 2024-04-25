function newstruct = preallocateStruct(inputstruct,set_size)
%PREALLOCATESTRUCT Summary of this function goes here
%   Detailed explanation goes here
if isempty(inputstruct)
    error('An existing structure must be provided to perform preallocation')
end
if numel(set_size)==1
    set_size = [1 set_size];
end
names = fieldnames(inputstruct);
command = @(name,size) sprintf('"%s",cell([%s])',name,num2str(size));
structcommand = cell(1,numel(names));
for i=1:numel(names)
    structcommand{i} = command(names{i},set_size);
end
structcommand = join(structcommand,',');
final = sprintf('struct(%s);',structcommand{1});
newstruct = eval(final);
end

