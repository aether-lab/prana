function [M1] = bootstrapping_dataremove(DSIZE,ENUM,MASK)
% --- Bootstrapping Data Removal ---

Nx = DSIZE(1);
Ny = DSIZE(2);
Nt = 1;

M1   = zeros(DSIZE);
RMAT = rand(Nx,Ny,Nt);
EN   = 0;

while sum(M1(:))/(Nx*Ny) < ENUM && EN < 1
    M1 = RMAT<EN;
    M1(MASK<0) = 0; 
    EN = EN + 0.005;    
end

M1 = double(M1);

end