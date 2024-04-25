function nccreatewrite(varargin)
%%Creates the NetCDF file directly, and these NetCDF files are readbale in GrADS.
%%Get rid of writing the long commands everytime to create the netcdf4 file
%%or to add the variable to the nc file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Ankur Kumar                       Monday; July 12, 2021
%                                              Version: 2.0
% https://www.nsstc.uah.edu/users/ankur.kumar/
% Doctoral Student
% The University of Alabama in Huntsville, USA
% Department of Atmospheric and Earth Science
% Graduate Research Assistant (GRA)
% National Aeronautics and Space Administration (NASA)
% Interagency Implementation and Advanced Concepts Team (IMPACT)
%
% Email: ankurk017@gmail.com
%        ankur.kumar@nsstc.uah.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function:
%           One can create a NetCDF file with the predefined MATLAB
%           function using nccreate and ncwrite. This function is a
%           combination of these two predefined functions. One can use this
%           to save the time and get rid of writing the same commands for
%           storing multi variables in a nc file.
%
% Syntax:
%           nccreatewrite('sample_file.nc','variable',{'lon','lat','lev'},randi(20,20,30,5))
%
% Inputs:
%           First input should be the name of the nc file in which you want
%           to store the data.
%           Second input should be the name of the variable in which you
%           want to store the specific variable.
%           Third argument should be in braces (not structure) which must
%           contains the number of variables as that of the size of the
%           data you want to store.
%           ex: If you want to store 'lat' whose dimensions is 5*1, then
%           third argument should be {'lat'}
%           The above is becasue of MATLAB also stores the variable
%           dimensions seprately, not in the variable list.
%           You can see this when you use ncdisp to see the listed
%           variables in nc file.
%           Fourth argument should be the data you want to write in nc
%           file.
%           Fifth argument is an optional argument, and if provided, it should 
%           be the attributes of the variable.
% Example:
% 
%         random_data = normrnd(3,10,[20,30,5])+273;
% 
%         filename='sample.nc';
% 
%         delete(filename)
% 
%         lon_att.standard_name='longitude';
%         lon_att.long_name='longitude';
%         lon_att.units='degrees_east';
%         lon_att.axis='X';
%         nccreatewrite(filename,'lon',{'lon'},[65:84],lon_att)
% 
% 
%         lat_att.standard_name='latitude';
%         lat_att.long_name='latitude';
%         lat_att.units='degrees_north';
%         lat_att.axis='Y';
%         nccreatewrite(filename,'lat',{'lat'},[1:30],lat_att)
% 
% 
%         lev_att.long_name='generic';
%         lev_att.units='level';
%         lev_att.axis='Z';
%         nccreatewrite(filename,'lev',{'lev'},[1:5],lev_att)
% 
% 
%         nccreatewrite(filename,'TC',{'lon','lat','lev'},random_data)
% 
%         ncdisp(filename)
%
% 
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               ankur.kumar@nsstc.uah.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>=6
    error('Only 5 argument needed. Type help nccreatewrite to get help of this function')
end
filename=varargin{1};
variable_name=varargin{2};
var=varargin{3};
data=varargin{4};
if nargin==4
    warning('You are not writing any Attributes to %s variable %s ',variable_name)
elseif nargin==5
    lon_att=varargin{5};
    att_fields=[fieldnames(lon_att) struct2cell(lon_att)];
end
data_size=size(data);
data_size_length=length(data_size);
if ~isempty(find(data_size==1))
    data_size_length=data_size_length-1;
    data_size=data_size(data_size~=1);
end
variable_length_dim=length(var);
if data_size_length~=variable_length_dim
    error('Dimensions of the data which you want to write is not matching with the given number of arguments for storing the dimensions of data.')
end
dimension_att=[];
for i=1:data_size_length
    dimension_att_temp={var{i},data_size(i)};
    dimension_att=[dimension_att dimension_att_temp];
end
nccreate(filename,variable_name,...
    'Dimensions',dimension_att,...
    'format','netcdf4')
ncwrite(filename,variable_name,data);
if nargin==5
    for att_size = 1: size(att_fields,1)
        ncwriteatt(filename,variable_name,   att_fields{att_size,1},   att_fields{att_size,2});
    end
end
% disp_var=sprintf('| Succesfully written %s in file %s |',variable_name,filename);
% disp('-------------------------------------------')
% disp(disp_var)
% disp('--------------------------------------------')
end