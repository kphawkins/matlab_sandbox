function [tk,xk] = min_jerk_knots(a, d, dx, min_dt)
a = a(end:-1:1);
x_i = polyval(a,0);
xk(1) = x_i;
tk(1) = 0;
k = 1;
while 1
    good_knots = [d, tk(k)+min_dt];
    for msign = [-1, 1]
        a0 = a(end)-(polyval(a,tk(k))+msign*dx);
        aptest = a;
        aptest(end) = a0
        roots(aptest)
        for r = roots(aptest)'
            if imag(r) == 0 && r > tk(k)
                good_knots(end+1) = r;
            end
        end
    end
    tk(k+1) = min(good_knots);
    good_knots
    k = k + 1;
    xk(k) = polyval(a,tk(k));
    if tk(k) == d
        break
    end
end
