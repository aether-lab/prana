function F=poly_3xy_123zmod(a,allx1data,allXdata,orderz)
% function F=poly_3xy_123z(a,data)

% Writen by M. Brady
% Edited and Commented by S. Raben

x1=allx1data(:,1);
x2=allx1data(:,2);
x3=allx1data(:,3);


if orderz==1                % cubic xy, linear z
    Fpoly=a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
        a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
        a(10).*x1.^3 + a(11).*x1.^2.*x2 + a(12).*x1.*x2.^2 +...
        a(13).*x2.^3 + a(14).*x1.^2.*x3 + a(15).*x1.*x2.*x3 +...
        a(16).*x2.^2.*x3;

elseif orderz==2            % cubic xy, quadratic z
    Fpoly=a(1) + a(2).*x1 + a(3).*x2 + a(4).*x3 + a(5).*x1.^2 +...
        a(6).*x1.*x2 + a(7).*x2.^2 + a(8).*x1.*x3 + a(9).*x2.*x3 +...
        a(10).*x3.^2 + a(11).*x1.^3 + a(12).*x1.^2.*x2 + a(13).*x1.*x2.^2 +...
        a(14).*x2.^3 + a(15).*x1.^2.*x3 + a(16).*x1.*x2.*x3 +...
        a(17).*x2.^2.*x3 + a(18).*x1.*x3.^2 + a(19).*x2.*x3.^2;    
    
else             % cubic xyz
    
    
end

F=allXdata-Fpoly;

