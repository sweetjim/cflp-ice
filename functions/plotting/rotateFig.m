function out = rotateFig(fig,dir)
pos = fig.Position;

corners = [pos(1) pos(2);           % Bottom left
    pos(1)+pos(3) pos(2);           % Top left
    pos(1) pos(2)+pos(4);           % Bottom right
    pos(1)+pos(3) pos(2)+pos(4)];   % Top right

screen = get(0,'ScreenSize');

corners(:,1) = corners(:,1)/screen(3);
corners(:,2) = corners(:,2)/screen(4);
switch dir
    case 'sideways'
        if corners(4,:)>.5          % If top right corner is in right half 
            
        else
            
        end
    case 'upright'
        
end

end

