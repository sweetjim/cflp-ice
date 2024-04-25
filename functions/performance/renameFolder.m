function renameFolder(root,oldname,newname)
oldf = fullfile(root,oldname);
if isfolder(oldf)
    newf = strrep(oldf,oldname,newname);
    mkdir(newf)
    movefile(fullfile(oldf,'*'),newf)
    rmdir(oldf)
end
end