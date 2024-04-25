function newROI = ROIrot90(ROI,image,rotation)
rotation = mod(rotation,4);
[Y,X]   = size(image);
Y       = Y-1; 
X       = X-1;
x       = ROI(1);
y       = ROI(2);
Dx      = ROI(3);
Dy      = ROI(4);

if ~rotation
    newROI = ROI;
    return
end

switch rotation
    case 1
        newROI = [y X-(x+Dx)+3 Dy Dx];
    case 2
        newROI = [X-(x+Dx)+3 Y-(y+Dy)+3 Dx Dy];
    case 3
        newROI = [Y-(y+Dy)+3 x Dy Dx];
end

end

