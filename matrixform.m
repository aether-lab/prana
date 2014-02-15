function [X,Y,U,V,Eval,C,D]=matrixform(x,y,u,v,eval,c,d)
% --- Vector to Matrix Subfunction ---

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
N=length(x);

%initialize matrices
U=nan(length(b),length(a),size(u,2));
V=nan(length(b),length(a),size(v,2));
Eval=-1*ones(length(b),length(a),size(eval,2));

%generate grid matrix
[X,Y]=meshgrid(a,b);

%generate variable matrices (nans where no data available)
for i=1:size(U,3)
    for n=1:N
        I=find(b==y(n));
        J=find(a==x(n));
        U(I,J,i) = u(n,i);
        V(I,J,i) = v(n,i);
        Eval(I,J,i) = eval(n);
    end
end
if ~isempty(c)
    C=nan(length(b),length(a),size(c,2));
    for i=1:size(c,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            C(I,J,i)=c(n,i);
        end
    end
else
    C=[];
end
if ~isempty(d)
    D=nan(length(b),length(a),size(d,2));
    for i=1:size(d,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            D(I,J,i)=d(n,i);
        end
    end
else
    D=[];
end

end