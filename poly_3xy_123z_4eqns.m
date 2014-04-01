
function F=poly_3xy_123z_4eqns(xyz,eqnsolvedata);
% this function will find x,y,z from the four equations
% from the known image coordinates (X1 Y1), (X2 Y2)
% xyz is the (1x3) vector of object (x y z) coordinates
% this can be used to map the surface point

X1=eqnsolvedata.X1Y1X2Y2(1);    
Y1=eqnsolvedata.X1Y1X2Y2(2);
X2=eqnsolvedata.X1Y1X2Y2(3);
Y2=eqnsolvedata.X1Y1X2Y2(4);
a1=eqnsolvedata.

orderz=eqnsolvedata.orderz;

if orderz==1                % cubic xy, linear z
    
    %F=
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
% disp(['a = ' num2str(a)]);
F=data.allXdata-Fpoly;

