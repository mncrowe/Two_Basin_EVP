function [H, Hx, Hy, Mask] = CreateBasins(x, y, n, b, sigma)
% Creates a n-basin geometry.
%
% Inputs:
% - x, y: grids, vectors of size (Nx, 1) and (1, Ny)
% - n: number of basins (default: 2)
% - b: basin steepness (default: 1)
% - sigma: list of sigma values [sigma_0 sigma_1 sigma_2 ...] (default:
%          [0.6 -0.5 -0.5])

arguments
    x (:,1) double
    y (1,:) double
    n (1,1) double     = 2
    b (1,1) double     = 1
    sigma (:,1) double = [0.6 -0.5 -0.5]
end

x(x >= x(end)/2) = x(x >= x(end)/2) - x(end);
y(y >= y(end)/2) = y(y >= y(end)/2) - y(end);

syms xs ys
assume(xs, "real")
assume(ys, "real")

z = xs + 1i*ys;

s_sym = 1/2*log(real((z^n - 1)*(z'^n-1)));
sx_sym = diff(s_sym, xs);
sy_sym = diff(s_sym, ys);

s_func = matlabFunction(s_sym);
sx_func = matlabFunction(sx_sym);
sy_func = matlabFunction(sy_sym);

s = s_func(x, y);
sx = sx_func(x, y);
sy = sy_func(x, y);

H = exp(-b*s);
Hp = -b*exp(-b*s);

theta = atan2(y, x);
domain = mod(floor((theta+pi/n)*n/(2*pi)), n);

for i = 1:n
    
    H((domain == (i-1)) & (s < sigma(i+1))) = exp(-b*sigma(i+1));
    Hp((domain == (i-1)) & (s < sigma(i+1))) = 0;

end

H(s > sigma(1)) = exp(-b*sigma(1));
Hp(s > sigma(1)) = 0;

Hx = sx .* Hp;
Hy = sy .* Hp;

Mask = 0*s;
Mask(s > sigma(1)) = 1;

end