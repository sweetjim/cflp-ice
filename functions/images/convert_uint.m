function image = convert_uint(image)
%CONVERT_UINT 
if max(image,[],'all')>256&&contains(class(image),'uint16')
    image = (image)*256/65535;
end
image = double(image);
end

