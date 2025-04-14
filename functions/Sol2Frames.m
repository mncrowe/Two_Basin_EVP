function Sol2Frames(psi, x, y, t, i, res, framename, N, s0)

    G = zeros(res(1), res(2), length(t));

    x_res = linspace(x(1), x(end), res(1));
    y_res = linspace(y(1), y(end), res(2));

    if nargin > 7

        syms xs ys
        assume(xs, "real")
        assume(ys, "real")

        z = xs + 1i*ys;

        s_sym = 1/2*log(real((z^N - 1)*(z'^N-1)));
        s_func = matlabFunction(s_sym);

        s = s_func(x_res', y_res);

        M = zeros(res);

        M(s > s0) = NaN;
        M(s <= s0) = 0;

    else

        M = 0;

    end

    for j = 1:length(t)
        F = real(psi(i, t(j)));
        G(:, :, j) = interp_grid(F, x, y, x_res, y_res) + M;
    end

    m = max(max(G,[],"all"), -min(G,[],"all"));

    save_frames(G, framename, 'png', 0, [-m m]);

end