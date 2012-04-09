function [locations] = circlecheck(locations,Icir)
[ymaxy,xmaxx] = size(Icir);
C_x = mean(locations(:,2));
C_y = mean(locations(:,1));
Cir = sqrt(((locations(:,1)-C_x).^2)+((locations(:,2)-C_y).^2));

while sum(abs(diff(Cir))) == 0
    if locations(1,1)-1 >= 1 && Icir(locations(1,1)-1,locations(1,2)) ~= 0
        locations(1,1) = locations(1,1)-1;
    elseif locations(2,1)+1 <= ymaxy  && Icir(locations(2,1)+1,locations(2,2)) ~= 0
        locations(2,1) = locations(2,1)+1;
    elseif locations(3,2)-1 >= 1 && Icir(locations(3,1),locations(3,2)-1) ~= 0
        locations(3,2) = locations(3,2)-1;
    elseif locations(4,2)+1 <= xmaxx && Icir(locations(4,1),locations(4,2)+1) ~= 0
        locations(4,2) = locations(4,2)+1;
    else
        %keyboard
    end
    
    C_x = mean(locations(:,2));
    C_y = mean(locations(:,1));
    Cir = sqrt(((locations(:,2)-C_x).^2)+((locations(:,1)-C_y).^2));
end
end
