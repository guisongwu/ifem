function ul = extend_mid(u)
    n = size(u, 1);
    ul = [u; u];
    ul(n+1:2*n, :) = .5 * (u(1:n, :) + [u(2:n, :); u(1, :)]);
end
