function F=campinhole_LSfun(a,data)
% F=campinhole_LSfun(a,data)
% solves the camera pinhole problem (linear), see Ng and Zhang, 2006, Exp. 
% Fluids., 34, pp. 484-493.  this solves equation 8 in that paper.  See 
% also Willert, exp. in fluids. (2006), n.41, pp.135., set the constraint 
% of a(12)=1, so the trivial solution is not solved for.
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

%  World target coord. for X,Y, and Z, x1,x2,x3 respectively.
x1=data.allxdata(:,1);
x2=data.allxdata(:,2);
x3=data.allxdata(:,3);

% Pixel target coordinates in X and Y direction, X1,X3 respectively.
X1=data.allXdata(:,1);
X2=data.allXdata(:,2);

% This equation comes from Willert 2006 eqn (2).  
eqns1=(a(1).*x1+a(2).*x2+a(3).*x3+a(4))-X1.*(a(9).*x1+a(10).*x2+a(11).*x3+1);
eqns2=(a(5).*x1+a(6).*x2+a(7).*x3+a(8))-X2.*(a(9).*x1+a(10).*x2+a(11).*x3+1);

% F=[eqns1;eqns2;eqns3]
F=[eqns1;eqns2];