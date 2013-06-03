function [] = test_time()
for i = 1:10000
    q = [0, 2*(rand(1,6)-0.5)*pi*2, pi];
    t = [0, cumsum(2*rand(1,9)+0.5)];
    qd_i = 0;
    qdd_i = 0;
    qd_f = 0;
    qdd_f = 0;


    h = t(2:end)-t(1:end-1);
    N = numel(t)-1;
    q = [q(1), 0, q(2:end-1), 0, q(end)]; % add in knot placeholders
    [a,b,c,d] = make_tri(h,q,qd_i,qdd_i,qd_f,qdd_f); % create tridiagonal matrix problem
    X = TDMAsolver(a,b,c,d); % solve tridiagonal matrix problem
end
