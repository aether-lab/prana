function F=campinholemod_LSfun(a,data)
% F=campinholemod_LSfun(a,data)
% solves the camera pinhole problem for a single component, see Ng and 
% Zhang, 2006, Exp. Fluids., 34, pp. 484-493.  this solves equation 8 
% in that paper.
% 
% This function is called by fitcameramodels.m
% 
% Inputs:
%    data:          Contains a file structure that holds the world and
%                   pixel corrdinates for the target markers.
%    data.allxdata: World coordinates for the target markers
%    data.allXdata: Pixel coordinates for the target markers
%    a:             Polynomial Coef.  There are only 11 coef. entered in
%                   because the 12th is set manually to 1.
% Output:
%    F:             Residual from the minimization on the pinhole camera 
%                   model equations. F(1) is camera 1, F(2) is camera 2.

% Writen by M. Brady
% Edited and Commented by S. Raben

x1=data.allxdata(:,1);
x2=data.allxdata(:,2);
x3=data.allxdata(:,3);

X=data.allXdata;


eqns1=a(1).*x1+a(2).*x2+a(3).*x3+a(4)-X;
% eqns2=a(5).*x1+a(6).*x2+a(7).*x3+a(8)-X2;
% eqns3=a(9).*x1+a(10).*x2+a(11).*x3+a(12)-1;

% eqns1=(X1*a(9)-a(1)).*x1+(X1*a(10)-a(2)).*x2+(X1*a(11)-a(3)).*x3 - a(4)+X1*a(12);
% eqns2=(X2*a(9)-a(5)).*x1+(X2*a(10)-a(6)).*x2+(X2*a(11)-a(7)).*x3 - a(8)+X2*a(12);

F=[eqns1];

% F=[eqns1;eqns2]