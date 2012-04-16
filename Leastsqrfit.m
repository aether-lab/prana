function [x_c,y_c,D,P,E,Meth] = Leastsqrfit(I_in,method_in,sigma_in)

%Desicription
%
%
%
%







%Options for the lsqnonlin solver
options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
    'Display','off','DiffMinChange',1e-7,'DiffMaxChange',1);%,'LevenbergMarquardt','off');%,'LargeScale','off');

% Find the center Max and all of the saturated points
[locxy_in(:,1) locxy_in(:,2)] = find(I_in == max(I_in(:)));
max_locxy_in(1) = round(median(locxy_in(:,1)));
max_locxy_in(2) = round(median(locxy_in(:,2)));

% If there are not enough points don't use the method
if nnz(I_in) - numel(find(I_in == max(I_in(:)))) + 1 < 5
    x_c = 1;
    y_c = 1;
    D   = NaN;
    P   = NaN;
    E   = NaN;
    Meth = 0;
%     fprintf('Not Enough Points for Least Squres Fit\n')
    return
end

%Removes negitive values and puts in zeros
Ils = I_in;
Ils(I_in<0) = 0;

[x_c,y_c,D,P,E,Meth] = Leastsqrmethods(Ils,method_in,sigma_in,options,locxy_in,max_locxy_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
