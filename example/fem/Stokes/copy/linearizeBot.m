function ul = linearizeBot(u)
    n = size(u, 1);
    assert(mod(n, 2) == 0);

    n2 = n/2;
    ul = u;
    ul(n2+1:n, :) = .5 * (u(1:n2, :) + [u(2:n2, :); u(1, :)]);
end


