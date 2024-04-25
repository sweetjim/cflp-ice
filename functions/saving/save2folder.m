function save2folder(path,new_struct,struct_name)
root        = fullfile(cd,'results');
folder      = fullfile(root,path);
eval(sprintf('%s=%s',struct_name,'new_struct'))
save(folder,struct_name,'-append')
end

