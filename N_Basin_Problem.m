% Solves the N-basin problem for Rossby wave modes

addpath("functions\misc\", "functions\")

% Define problem parameters:

N = 2;          % number of basins
n = 10;         % number of modes to find
w0 = 0.25;      % frequency guess

sigma = [0.6 -0.5*ones(1, N)];    % sigma values on boundary contours

BC = 1;         % boundary condition type: 1 - psi = 0 on s = sigma_0
                %                          0 - periodic in x, y

% Define grid parameters:

Nx = 121;
Ny = 121;

Lx = 4.0;
Ly = 2.0;

% Define output parameters:

save_movies = 1;    % if set to 1, output is saved as movie files
show_mask = 1;      % if set to 1, region sigma > sigma_0 shown in grey
N_res = [1000 500];  % number of pixels in each direction in output movies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       End of User Input       %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define helper functions:

to_vec = @(F) reshape(F, [Nx*Ny 1]);
to_mat = @(F) reshape(F, [Nx Ny]);

% Create grids:

[Mx, x] = grid_spectral(2, Nx, [0 Lx]);
[My, y] = grid_spectral(2, Ny, [0 Ly]);

% Construct topography profile:

[h, hx, hy, Mask] = CreateBasins(x, y', N, 1, sigma);

% Construct eigenvalue problem:

disp('Building EVP...')

Ix = speye(Nx);
Iy = speye(Ny);

Hx = sparse(1:Nx*Ny, 1:Nx*Ny, to_vec(hx ./ h));
Hy = sparse(1:Nx*Ny, 1:Nx*Ny, to_vec(hy ./ h));

Dx = kron(Iy, Mx);
Dy = kron(My, Ix);

D2x = kron(Iy, Mx^2);
D2y = kron(My^2, Ix);

A = D2x + D2y - Hx * Dx - Hy * Dy;
B = Hx * Dy - Hy * Dx;

% Apply BCs, psi and normal derivatives are continuous across boundaries:

% BC at x = -Lx/2:
A(1:Nx:Nx*Ny, :) = 0;
B(1:Nx:Nx*Ny, :) = kron(Iy, [1 zeros(1, Nx-2) -1]);

% BC at x = Lx/2:
A(Nx:Nx:Nx*Ny, :) = 0;
B(Nx:Nx:Nx*Ny, :) = kron(Iy, Mx(1, :) - Mx(Nx, :));

% BC at y = -Ly/2:
A(1:Nx, :) = 0;
B(1:Nx, :) = kron([1 zeros(1, Ny-2) -1], Ix);

% BC at y = Ly/2:
A(Nx*(Ny-1)+1:Nx*Ny, :) = 0;
B(Nx*(Ny-1)+1:Nx*Ny, :) = kron(My(1, :) - My(Ny, :), Ix);

% If BC = 0, apply gauge condition (psi = 0) at domain centre (x,y = 0):

if BC == 0
    Nc = (Nx*Ny + 1) / 2;
    A(Nc, :) = 0;
    B(Nc, :) = [zeros(1, Nc-1) 1 zeros(1, Nc-1)];
end

% if BC = 1, set all points outside s = sigma_0 to zero:

if BC == 1
    D = sparse(diag(1:Nx*Ny));
    Mask_i = (to_vec(Mask)==1);
    A(Mask_i, :) = 0;
    B(Mask_i, :) = D(Mask_i, :);
end

%row_index = reshape(1:Nx*Ny, [Nx Ny]); % gives array showing which points
%in x,y grid correspond to which rows in A and B.

% solve for s = i\omega and \hat{\psi}:

disp('Solving EVP...')

[s, evec] = Gen_EVP(A, B, n, w0*1i);

w = real(s/1i);

[w, I] = sort(w, 'descend');
evec = evec(:, I);

psih = zeros(Nx, Ny, n);

for i = 1:n
    psih(:, :, i) = DomainShift(to_mat(evec(:,i)),x, y);
end

[h, x, y] = DomainShift(h, x, y);

% Costruct psi = Re[\hat{\psi} exp(-i\omega t)]:

psi = @(i, t) cos(w(i)*t) * real(psih(:,:,i)) + sin(w(i)*t) * imag(psih(:,:,i));

% Save output as a movie for each mode:

if save_movies == 1

    if show_mask == 1
        Mask = DomainShift(Mask, x, y);
    else
        Mask = zeros(Nx, Ny);
    end

    % create output directories:

    if ~exist("frames", "dir"); mkdir("frames"); end
    if ~exist("mov", "dir"); mkdir("mov"); end
    
    % Clear output directories:
    
    system(['del mov\mov_' num2str(N) '_basin* /Q 2>nul']);
    system('del frames\* /Q 2>nul');
    
    % Create frames and movie files:
    
    for i = 1:n
        disp(['Making movie ' num2str(i) ' of ' num2str(n) '...'])
        Sol2Frames(psi, x, y, linspace(0,4*pi/w(i),101), i, N_res, ['frames/frame_' num2str(i)], N, sigma(1));
        mov_name = ['mov/mov_' num2str(N) '_basin_mode_' num2str(i) '_w_' num2str(round(w(i), 3)) '.mp4'];
        system(['ffmpeg.exe -framerate 20 -i frames/frame_' num2str(i) '_%03d.png ' mov_name ' >nul 2>nul']);
    end
    
    close all

end
