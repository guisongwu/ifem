%% functions

syms x y y1 slope;


y1 = y + x * slope;

u = (1-y1)^2;

fprintf('Grad\n');
diff(u, x)
diff(u, y)

fprintf('-Lap\n');
-diff(u, x, 2) - diff(u, y, 2)



uu = 3/4 * (1-y) - 1/2*(1-y)^2;
subs(-diff(uu, y) + uu, y, 0)
