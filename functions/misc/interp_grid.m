function G = interp_grid(F,x1,y1,x2,y2,method,c)
% Interpolates the array F from the grid (x1,y1) to the grid (x2,y2)
% through the first two dimensions.
% - F: input field, n dimensional array for n >= 2
% - (x1,y1): input grids, must match size of first 2 dimensions of F
% - (x2,y2): output grids, F is interpolated onto these
% - method: interpolation method;
%       1 - quadratic interpolation (default)
%       2 - uses interp2_xy.m
%       3 - uses interp2 with 'spline'
% - c: optional arguments for interp2_xy, enter as [n m]
%
% Note: the quadratic interpolation assumes that the function can be
% modelled as quadratic in both the x and y directions using the nine (x,y)
% points closest to the interpolation point. The coefficients of the fit
% are found explicitly rather than just inverting a linear system as this
% method can be vectorised easily and is faster for large grids.
% The interpolation is:
%   F = N11 + N21*x + N12*y + N22*xy + N31*x^2 + N13*y^2 ... 
%                                   + N23*xy^2 + N32*x^2y + N33*x^2y^2

% set undefined parameters
if nargin < 6; method = 1; end
if nargin < 7; c = [1 2]; end

% convert grids to column vectors
x1 = reshape(x1,[length(x1) 1]); y1 = reshape(y1,[length(y1) 1]);
x2 = reshape(x2,[length(x2) 1]); y2 = reshape(y2,[length(y2) 1]);

% check for consistency and reshape to a 3D array
s = size(F);
if length(x1) ~= s(1) || length(y1) ~= s(2); error('Size of input inconsistent with grid'); end
if length(s) > 2; F = reshape(F,[s(1) s(2) prod(s(3:end))]); end

% check that (x2,y2) lies within the ranges of (x1,y1)
if min(x2) < min(x1) || max(x2) > max(x1) || min(y2) < min(y1) || max(y2) > max(y1); error('interpolation points outside input domain'); end

% define output array
G = zeros(length(x2),length(y2),prod(s(3:end)));

% find closest point array (method 1 only)
if method == 1
    Nx = length(x1); Ny = length(y1);
    [~,ix] = min(abs(x1-x2')); ix = max(2,min(ix,Nx-1))';
    [~,iy] = min(abs(y1-y2')); iy = max(2,min(iy,Ny-1))';
    x_1 = x1(ix-1); x_2 = x1(ix); x_3 = x1(ix+1);
    y_1 = y1(iy-1)'; y_2 = y1(iy)'; y_3 = y1(iy+1)';

    Nx11 = -x_2.*x_3./(x_1-x_2)./(x_3-x_1); Nx12 = -x_1.*x_3./(x_1-x_2)./(x_2-x_3); Nx13 = -x_1.*x_2./(x_3-x_1)./(x_2-x_3);
    Nx21 = (x_2+x_3)./(x_1-x_2)./(x_3-x_1); Nx22 = (x_1+x_3)./(x_1-x_2)./(x_2-x_3); Nx23 = (x_1+x_2)./(x_3-x_1)./(x_2-x_3);
    Nx31 = -1./(x_1-x_2)./(x_3-x_1); Nx32 = -1./(x_1-x_2)./(x_2-x_3); Nx33 = -1./(x_3-x_1)./(x_2-x_3);
    Ny11 = -y_2.*y_3./(y_1-y_2)./(y_3-y_1); Ny12 = -y_1.*y_3./(y_1-y_2)./(y_2-y_3); Ny13 = -y_1.*y_2./(y_3-y_1)./(y_2-y_3);
    Ny21 = (y_2+y_3)./(y_1-y_2)./(y_3-y_1); Ny22 = (y_1+y_3)./(y_1-y_2)./(y_2-y_3); Ny23 = (y_1+y_2)./(y_3-y_1)./(y_2-y_3);
    Ny31 = -1./(y_1-y_2)./(y_3-y_1); Ny32 = -1./(y_1-y_2)./(y_2-y_3); Ny33 = -1./(y_3-y_1)./(y_2-y_3);
    
    N11 = Nx11.*Ny11 + x2.*Nx21.*Ny11 + Nx11.*Ny21.*y2' + x2.^2.*Nx31.*Ny11 + Nx11.*Ny31.*y2'.^2 + x2.*Nx21.*Ny21.*y2' + x2.^2.*Nx31.*Ny21.*y2' + x2.*Nx21.*Ny31.*y2'.^2 + x2.^2.*Nx31.*Ny31.*y2'.^2;
    N12 = Nx11.*Ny12 + x2.*Nx21.*Ny12 + Nx11.*Ny22.*y2' + x2.^2.*Nx31.*Ny12 + Nx11.*Ny32.*y2'.^2 + x2.*Nx21.*Ny22.*y2' + x2.^2.*Nx31.*Ny22.*y2' + x2.*Nx21.*Ny32.*y2'.^2 + x2.^2.*Nx31.*Ny32.*y2'.^2;
    N21 = Nx12.*Ny11 + x2.*Nx22.*Ny11 + Nx12.*Ny21.*y2' + x2.^2.*Nx32.*Ny11 + Nx12.*Ny31.*y2'.^2 + x2.*Nx22.*Ny21.*y2' + x2.^2.*Nx32.*Ny21.*y2' + x2.*Nx22.*Ny31.*y2'.^2 + x2.^2.*Nx32.*Ny31.*y2'.^2;
    N22 = Nx12.*Ny12 + x2.*Nx22.*Ny12 + Nx12.*Ny22.*y2' + x2.^2.*Nx32.*Ny12 + Nx12.*Ny32.*y2'.^2 + x2.*Nx22.*Ny22.*y2' + x2.^2.*Nx32.*Ny22.*y2' + x2.*Nx22.*Ny32.*y2'.^2 + x2.^2.*Nx32.*Ny32.*y2'.^2;
    N31 = Nx13.*Ny11 + x2.*Nx23.*Ny11 + Nx13.*Ny21.*y2' + x2.^2.*Nx33.*Ny11 + Nx13.*Ny31.*y2'.^2 + x2.*Nx23.*Ny21.*y2' + x2.^2.*Nx33.*Ny21.*y2' + x2.*Nx23.*Ny31.*y2'.^2 + x2.^2.*Nx33.*Ny31.*y2'.^2;
    N13 = Nx11.*Ny13 + x2.*Nx21.*Ny13 + Nx11.*Ny23.*y2' + x2.^2.*Nx31.*Ny13 + Nx11.*Ny33.*y2'.^2 + x2.*Nx21.*Ny23.*y2' + x2.^2.*Nx31.*Ny23.*y2' + x2.*Nx21.*Ny33.*y2'.^2 + x2.^2.*Nx31.*Ny33.*y2'.^2;
    N33 = Nx13.*Ny13 + x2.*Nx23.*Ny13 + Nx13.*Ny23.*y2' + x2.^2.*Nx33.*Ny13 + Nx13.*Ny33.*y2'.^2 + x2.*Nx23.*Ny23.*y2' + x2.^2.*Nx33.*Ny23.*y2' + x2.*Nx23.*Ny33.*y2'.^2 + x2.^2.*Nx33.*Ny33.*y2'.^2;
    N23 = Nx12.*Ny13 + x2.*Nx22.*Ny13 + Nx12.*Ny23.*y2' + x2.^2.*Nx32.*Ny13 + Nx12.*Ny33.*y2'.^2 + x2.*Nx22.*Ny23.*y2' + x2.^2.*Nx32.*Ny23.*y2' + x2.*Nx22.*Ny33.*y2'.^2 + x2.^2.*Nx32.*Ny33.*y2'.^2;
    N32 = Nx13.*Ny12 + x2.*Nx23.*Ny12 + Nx13.*Ny22.*y2' + x2.^2.*Nx33.*Ny12 + Nx13.*Ny32.*y2'.^2 + x2.*Nx23.*Ny22.*y2' + x2.^2.*Nx33.*Ny22.*y2' + x2.*Nx23.*Ny32.*y2'.^2 + x2.^2.*Nx33.*Ny32.*y2'.^2;
end

if method == 1
    G = F(ix-1,iy-1,:).*N11 + F(ix,iy-1,:).*N21 + F(ix-1,iy,:).*N12 + F(ix,iy,:).*N22 + F(ix+1,iy-1).*N31 + F(ix-1,iy+1,:).*N13 + F(ix+1,iy+1,:).*N33 + F(ix+1,iy,:).*N32 + F(ix,iy+1,:).*N23;
end
if method == 2
    for i = 1:prod(s(3:end))
        for ix = 1:length(x2)
            for iy = 1:length(y2)
                G(ix,iy,i) = interp2_xy(squeeze(F(:,:,i)),x1,y1,x2(ix),y2(iy),c(1),c(2));
            end
        end
    end
end
if method == 3
    for i = 1:prod(s(3:end))
        G(:,:,i) = interp2(x1,y1',F(:,:,i),x2,y2','spline');
    end
end

% reshape back to original number of dimensions
if length(s) > 2; G = reshape(G,[length(x2) length(y2) s(3:end)]); end

end

