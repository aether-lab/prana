function [W]=energyfilt(Nx,Ny,d,q)
% --- RPC Spectral Filter Subfunction ---

%     This file is part of prana, an open-source GUI-driven program for
%     calculating velocity fields using PIV or PTV.
%     Copyright (C) 2012  Virginia Polytechnic Institute and State
%     University
% 
%     prana is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


if numel(d) == 1
    d(2) = d;
end

%assume no aliasing
if nargin<4
    q = 0;
end

%initialize indices
[k1,k2]=meshgrid(-pi:2*pi/Ny:pi-2*pi/Ny,-pi:2*pi/Nx:pi-2*pi/Nx);

%particle-image spectrum
% Ep = (pi*255*d^2/8)^2*exp(-d^2*k1.^2/16).*exp(-d^2*k2.^2/16);
Ep = (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*k1.^2/16).*exp(-d(1)^2*k2.^2/16);

%aliased particle-image spectrum
% Ea = (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2+2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1+0*pi).^2/16).*exp(-d^2*(k2-2*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1+2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16)+...
%      (pi*255*d^2/8)^2*exp(-d^2*(k1-2*pi).^2/16).*exp(-d^2*(k2+0*pi).^2/16);
Ea = (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+0*pi).^2/16).*exp(-d(1)^2*(k2+2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+0*pi).^2/16).*exp(-d(1)^2*(k2-2*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1+2*pi).^2/16).*exp(-d(1)^2*(k2+0*pi).^2/16)+...
     (pi*255*(d(1)*d(2))/8)^2*exp(-d(2)^2*(k1-2*pi).^2/16).*exp(-d(1)^2*(k2+0*pi).^2/16);

%noise spectrum
En = pi/4*Nx*Ny;

%DPIV SNR spectral filter
W  = Ep./((1-q)*En+(q)*Ea);
W  = W'/max(max(W));


end
