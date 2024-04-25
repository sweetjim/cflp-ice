function files = readfolder
root    = uigetdir;
items   = dir(root);
files   = items(3:end);
end

