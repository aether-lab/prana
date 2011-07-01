function [W]=windowmask(N,R)
% --- Gaussian Window Mask Subfunction ---

% %generic indices
x  = -1:2/(N(1)-1):1;
y  = (-1:2/(N(2)-1):1)';
% 
% %gaussian window sizes
% px = (1.224*N(1)/R(1))^1.0172;
% py = (1.224*N(2)/R(2))^1.0172;
[px]=findwidth(R(1)/N(1));
[py]=findwidth(R(2)/N(2));
% 
% %generate 2D window
wx=exp(-px^2.*x.^2/2);
wy=exp(-py^2.*y.^2/2);

W  = wy*wx;