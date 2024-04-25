function ncfileProtocol(savepath,varname,varval)
%NCFILEPROTOCOL is a function that enables predefined Dimension and
%Attribute handling for creating and writing to NetCDF files.
%_________________________________________________________________
%_________________________________________________________________
% Input arguments:
%   'varname'       -   [string]
%                       Variable name.
%   'varval'        -   [nxmxk matrix, double]
%                       Variable matrix.
%_________________________________________________________________
%_________________________________________________________________
%%
[xs,zs]  = size(varval);
nccreate(savepath,varname,'Dimensions',{'xx',xs,'yy',zs,'time',Inf});
ncwrite(savepath,varname,varval);    
end

