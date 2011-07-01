function [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes,varargin)
% gridfit: estimates a surface on a 2d grid, based on scattered data
%          Replicates are allowed. All methods extrapolate to the grid
%          boundaries. Gridfit uses a modified ridge estimator to
%          generate the surface, where the bias is toward smoothness.
%
%          Gridfit is not an interpolant. Its goal is a smooth surface
%          that approximates your data, but allows you to control the
%          amount of smoothing.
%
% usage #1: zgrid = gridfit(x,y,z,xnodes,ynodes);
% usage #2: [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes);
% usage #3: zgrid = gridfit(x,y,z,xnodes,ynodes,prop,val,prop,val,...);
%
% Arguments: (input)
%  x,y,z - vectors of equal lengths, containing arbitrary scattered data
%          The only constraint on x and y is they cannot ALL fall on a
%          single line in the x-y plane. Replicate points will be treated
%          in a least squares sense.
%
%          ANY points containing a NaN are ignored in the estimation
%
%  xnodes - vector defining the nodes in the grid in the independent
%          variable (x). xnodes need not be equally spaced. xnodes
%          must completely span the data. If they do not, then the
%          'extend' property is applied, adjusting the first and last
%          nodes to be extended as necessary. See below for a complete
%          description of the 'extend' property.
%
%          If xnodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%  ynodes - vector defining the nodes in the grid in the independent
%          variable (y). ynodes need not be equally spaced.
%
%          If ynodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.

% set defaults
params.smoothness = 1;
params.interp = 'triangle';
params.regularizer = 'gradient';
params.solver = 'normal';
params.maxiter = [];
params.extend = 'warning';

% and check for any overrides
params = parse_pv_pairs(params,varargin);

% check the parameters for acceptability
% smoothness == 1 by default
if isempty(params.smoothness)
  params.smoothness = 1;
else
  if (params.smoothness<=0)
    error 'Smoothness must be real, finite, and positive.'
  end
end
% regularizer  - must be one of 4 options - the second and
% third are actually synonyms.
valid = {'springs', 'diffusion', 'laplacian', 'gradient'};
if isempty(params.regularizer)
  params.regularizer = 'diffusion';
end
ind = strmatch(lower(params.regularizer),valid);
if (length(ind)==1)
  params.regularizer = valid{ind};
else
  error(['Invalid regularization method: ',params.regularizer])
end

% interp must be one of:
%    'bilinear', 'nearest', or 'triangle'
% but accept any shortening thereof.
valid = {'bilinear', 'nearest', 'triangle'};
if isempty(params.interp)
  params.interp = 'triangle';
end
ind = strmatch(lower(params.interp),valid);
if (length(ind)==1)
  params.interp = valid{ind};
else
  error(['Invalid interpolation method: ',params.interp])
end

% solver must be one of:
%    'backslash', '\', 'symmlq', 'lsqr', or 'normal'
% but accept any shortening thereof.
valid = {'backslash', '\', 'symmlq', 'lsqr', 'normal'};
if isempty(params.solver)
  params.solver = '\';
end
ind = strmatch(lower(params.solver),valid);
if (length(ind)==1)
  params.solver = valid{ind};
else
  error(['Invalid solver option: ',params.solver])
end

% extend must be one of:
%    'never', 'warning', 'always'
% but accept any shortening thereof.
valid = {'never', 'warning', 'always'};
if isempty(params.extend)
  params.extend = 'warning';
end
ind = strmatch(lower(params.extend),valid);
if (length(ind)==1)
  params.extend = valid{ind};
else
  error(['Invalid extend option: ',params.extend])
end

% ensure all of x,y,z,xnodes,ynodes are column vectors,
% also drop any NaN data
x=x(:);
y=y(:);
z=z(:);
k = isnan(x) | isnan(y) | isnan(z);
if any(k)
  x(k)=[];
  y(k)=[];
  z(k)=[];
end
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% did they supply a scalar for the nodes?
if length(xnodes)==1
  xnodes = linspace(xmin,xmax,xnodes)';
  xnodes(end) = xmax; % make sure it hits the max
end
if length(ynodes)==1
  ynodes = linspace(ymin,ymax,ynodes)';
  ynodes(end) = ymax; % make sure it hits the max
end

xnodes=xnodes(:);
ynodes=ynodes(:);
dx = diff(xnodes);
dy = diff(ynodes);
nx = length(xnodes);
ny = length(ynodes);
ngrid = nx*ny;

% default for maxiter?
if isempty(params.maxiter)
  params.maxiter = min(10000,nx*ny);
end

% check lengths of the data
n = length(x);
if (length(y)~=n)||(length(z)~=n)
  error 'Data vectors are incompatible in size.'
end
if n<3
  error 'Insufficient data for surface estimation.'
end

% verify the nodes are distinct
if any(diff(xnodes)<=0)||any(diff(ynodes)<=0)
  error 'xnodes and ynodes must be monotone increasing'
end

% do we need to tweak the first or last node in x or y?
if xmin<xnodes(1)
    xnodes(1) = xmin;
end
if xmax>xnodes(end)
    xnodes(end) = xmax;
end
if ymin<ynodes(1)
    ynodes(1) = ymin;
end
if ymax>ynodes(end)
    ynodes(end) = ymax;
end

% only generate xgrid and ygrid if requested.
if nargout>1
  [xgrid,ygrid]=meshgrid(xnodes,ynodes);
end

% determine which cell in the array each point lies in
[junk,indx] = histc(x,xnodes);
[junk,indy] = histc(y,ynodes);
% any point falling at the last node is taken to be
% inside the last cell in x or y.
k=(indx==nx);
indx(k)=indx(k)-1;
k=(indy==ny);
indy(k)=indy(k)-1;

% interpolation equations for each point
tx = min(1,max(0,(x - xnodes(indx))./dx(indx)));
ty = min(1,max(0,(y - ynodes(indy))./dy(indy)));
ind = indy + ny*(indx-1);
% Future enhancement: add cubic interpolant
switch params.interp
  case 'triangle'
    % linear interpolation inside each triangle
    k = (tx > ty);
    L = ones(n,1);
    L(k) = ny;
    
    t1 = min(tx,ty);
    t2 = max(tx,ty);
    A = sparse(repmat((1:n)',1,3),[ind,ind+ny+1,ind+L], ...
       [1-t2,t1,t2-t1],n,ngrid);
    
  case 'nearest'
    % nearest neighbor interpolation in a cell
    k = round(1-ty) + round(1-tx)*ny;
    A = sparse((1:n)',ind+k,ones(n,1),n,ngrid);
    
  case 'bilinear'
    % bilinear interpolation in a cell
    A = sparse(repmat((1:n)',1,4),[ind,ind+1,ind+ny,ind+ny+1], ...
       [(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], ...
       n,ngrid);
    
end
rhs = z;

% Build regularizer. Add del^4 regularizer one day.
switch params.regularizer
  case 'springs'
    % zero "rest length" springs
    [i,j] = meshgrid(1:nx,1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = nx*(ny-1);
    stiffness = 1./dy;
    Areg = sparse(repmat((1:m)',1,2),[ind,ind+1], ...
       stiffness(j(:))*[-1 1],m,ngrid);
    
    [i,j] = meshgrid(1:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*ny;
    stiffness = 1./dx;
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny], ...
       stiffness(i(:))*[-1 1],m,ngrid)];
    
    [i,j] = meshgrid(1:(nx-1),1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*(ny-1);
    stiffness = 1./sqrt(dx(i(:)).^2 + dy(j(:)).^2);
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny+1], ...
       stiffness*[-1 1],m,ngrid)];
    
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind+1,ind+ny], ...
       stiffness*[-1 1],m,ngrid)];
    
  case {'diffusion' 'laplacian'}
    % thermal diffusion using Laplacian (del^2)
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));
    
    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
       -2./(dy2.*(dy1+dy2))],ngrid,ngrid);
    
    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));
    
    Areg = Areg + sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
       -2./(dx2.*(dx1+dx2))],ngrid,ngrid);
    
  case 'gradient'
    % Subtly different from the Laplacian. A point for future
    % enhancement is to do it better for the triangle interpolation
    % case.
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
      -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
      -2./(dx2.*(dx1+dx2))],ngrid,ngrid)];

end
nreg = size(Areg,1);

% Append the regularizer to the interpolation equations,
% scaling the problem first. Use the 1-norm for speed.
NA = norm(A,1);
NR = norm(Areg,1);
A = [A;Areg*(params.smoothness*NA/NR)];
rhs = [rhs;zeros(nreg,1)];

% solve the full system, with regularizer attached
switch params.solver
  case {'\' 'backslash'}
    % permute for minimum fill in for R (in the QR)
    p = colamd(A);
    zgrid=zeros(ny,nx);
    zgrid(p) = A(:,p)\rhs;
    
  case 'normal'
    % The normal equations, solved with \. Can be fast
    % for huge numbers of data points.
    
    % Permute for minimum fill-in for \ (in chol)
    APA = A'*A;
    p = symamd(APA);
    zgrid=zeros(ny,nx);
    zgrid(p) = APA(p,p)\(A(:,p)'*rhs);
    
  case 'symmlq'
    % iterative solver - symmlq - requires a symmetric matrix,
    % so use it to solve the normal equations. No preconditioner.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = symmlq(A'*A,A'*rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % SYMMLQ iterated MAXIT times but did not converge.
        warning(['Symmlq performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % SYMMLQ stagnated, successive iterates were the same
        warning 'Symmlq stagnated without apparent convergence.'
      otherwise
        warning(['One of the scalar quantities calculated in',...
          ' symmlq was too small or too large to continue computing.'])
    end
    
  case 'lsqr'
    % iterative solver - lsqr. No preconditioner here.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = lsqr(A,rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % lsqr iterated MAXIT times but did not converge.
        warning(['Lsqr performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % lsqr stagnated, successive iterates were the same
        warning 'Lsqr stagnated without apparent convergence.'
      case 4
        warning(['One of the scalar quantities calculated in',...
          ' LSQR was too small or too large to continue computing.'])
    end
    
end