function [u,v,eval,c]=vectorform(x,y,U,V,Eval,C)
% --- Matrix to Vector Subfunction ---
x=x(:);y=y(:);
%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x(:));

%initialize vectors
S=size(x(:));
u    = zeros(S);
v    = zeros(S);
eval = zeros(S);
if ~isempty(C)
    c = zeros(S);
else
    c = [];
end

%generate data vectors where data is available
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    u(n)    = U(I,J);
    v(n)    = V(I,J);
    eval(n) = Eval(I,J);
    if ~isempty(C)
        c(n)    = C(I,J);
    end
end