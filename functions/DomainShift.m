function [F, x, y] = DomainShift(F, x, y)
% Shifts a field and the corresponding grid to a domain centred at zero
%
% Note: x and y should have an ODD number of points

Nx = length(x);
Ny = length(y);

Nx2 = (Nx+1) / 2;
Ny2 = (Ny+1) / 2;

F = F([Nx2:Nx-1 1:Nx2], [Ny2:Ny-1 1:Ny2]);

x = [x(Nx2:Nx-1)-x(end); x(1:Nx2)];
y = [y(Ny2:Ny-1)-y(end); y(1:Ny2)];

end