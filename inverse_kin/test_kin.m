q = [ 4.07545758,  5.71643082, -4.57552159, -2.79061482, -3.17069678, 1.42865389];
arm_params = [0.1273, -0.612, -0.5723, 0.163941, 0.1157, 0.0922];
all_knowns = [x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, d1, a2, a3, d4, d5, d6];
arm_knowns = [d1, a2, a3, d4, d5, d6];
joints = [q1, q2, q3, q4, q5, q6];
Btests{1} = subs(Bmats(1), joints(1), q(1));
for i = 2:6
    Btests{i} = Btests{i-1} * subs(Bmats(i), joints(i), q(i));
    Btests{i} = subs(Btests{i}, arm_knowns, arm_params);
end
Btest = Btests{6}
Btest_act = vpa(Bb0*Btest*B6e)
mat_params = [Btest(1,1), Btest(2,1), Btest(3,1), Btest(1,2), Btest(2,2), Btest(3,2), Btest(1,3), Btest(2,3), Btest(3,3), Btest(1,4), Btest(2,4), Btest(3,4)];
known_vals = [mat_params, arm_params];
all_knowns2 = [all_knowns, q1, q2, q3, q4, q5, q6];
known_vals2 = [known_vals, q];
for i = 1:6
    qg(i) = vpa(subs(alleqssimp(i), all_knowns, known_vals))
    vpa(subs(alleqssimp(i), all_knowns2, known_vals2))
end
