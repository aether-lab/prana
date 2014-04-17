function [u,v,eval,c,d]=vectorform(x,y,U,V,Eval,C,D)
% --- Matrix to Vector Subfunction ---
x=x(:);y=y(:);
%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x(:));

imClass = 'double';

%initialize vectors
S=size(x(:));
u    = zeros(S,imClass);
v    = zeros(S,imClass);
eval = zeros(S,imClass);
if ~isempty(C)
    c = zeros(S,imClass);
    d = zeros(S,imClass);
else
    c = [];
    d = [];
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
        d(n)    = D(I,J);
    end
end

end
