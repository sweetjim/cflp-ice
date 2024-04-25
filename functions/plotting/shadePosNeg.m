function shadePosNeg(x1d,y1d,direction)
%SHADEPOSNEG 
if nargin<3
    direction = 'x';
end
%%
yneg = y1d;
ypos = y1d;

yneg(yneg>=0)=0;
ypos(ypos<=0)=0;
yneg([1 end])=0;
ypos([1 end])=0;
switch direction
    case 'x'
        patch(x1d,yneg,'r','FaceAlpha',.1,'EdgeAlpha',0)
        patch(x1d,ypos,'b','FaceAlpha',.1,'EdgeAlpha',0)
    case 'y'
        patch(yneg,x1d,'r','FaceAlpha',.1,'EdgeAlpha',0)
        patch(ypos,x1d,'b','FaceAlpha',.1,'EdgeAlpha',0)
end
end

