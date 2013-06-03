syms q1 q2 q3 q4 q5 q6 d1 a2 a3 d4 d5 d6
syms x1 x2 x3 y1 y2 y3 z1 z2 z3 p1 p2 p3
knowns = [x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, d1, a2, a3, d4, d5, d6];
known_vals = sym([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]);
if 0
B01 = [cos(q1), -sin(q1), 0, 0;
       sin(q1), cos(q1), 0, 0;
       0, 0, 1, d1;
       0, 0, 0, 1];
B12 = [cos(q2), -sin(q2), 0, 0;
       0, 0, -1, 0;
       sin(q2), cos(q2), 0, 0;
       0, 0, 0, 1];
B23 = [cos(q3), -sin(q3), 0, a2;
       sin(q3), cos(q3), 0, 0;
       0, 0, 1, 0;
       0, 0, 0, 1];
B34 = [cos(q4), -sin(q4), 0, a3;
       sin(q4), cos(q4), 0, 0;
       0, 0, 1, d4;
       0, 0, 0, 1];
B45 = [cos(q5), -sin(q5), 0, 0;
       0, 0, -1, -d5;
       sin(q5), cos(q5), 0, 0;
       0, 0, 0, 1];
B56 = [cos(q6), -sin(q6), 0, 0;
       0, 0, 1, d6;
       sin(q6), cos(q6), 0, 0;
       0, 0, 0, 1];
Bmats = {B01, B12, B23, B34, B45, B56};
end
if 0
    Bb0 = sym([-1, 0, 0, 0; 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]);
    B6e = sym([0, -1, 0, 0; 0, 0, -1, 0; 1, 0, 0, 0; 0, 0, 0, 1]);
    qsym = [q1, q2, q3, q4, q5, q6];
    asym = [0, a2, a3, 0, 0, 0];
    dsym = [d1, 0, 0, d4, d5, d6];
    lsym = sym([pi/2, 0, 0, pi/2, -pi/2, 0]);
    Bmats = {};
    for i = 1:6
        qi = qsym(i);
        ai = asym(i);
        di = dsym(i);
        li = lsym(i);
        B = [cos(qi), -sin(qi)*cos(li),  sin(qi)*sin(li), ai*cos(qi);
             sin(qi),  cos(qi)*cos(li), -cos(qi)*sin(li), ai*sin(qi);
                   0,          sin(li),          cos(li),         di;
                   0,                0,                0,          1]
        Bmats{i} = B;
    end
end
B06 = [x1, y1, z1, p1;
       x2, y2, z2, p2;
       x3, y3, z3, p3;
       0, 0, 0, 1];


if 0

numvars = {};
Bsols = {};
for i = 1:6
    i
    Bleft = eye(4);
    for j = 1:i-1
        Bleft = Bmats{j}^-1 * Bleft;
    end
    Bright = eye(4);
    for j = i:6
        Bright = Bright * Bmats{j}^-1;
    end
    Bfull = Bleft * B06 - Bright;
    Bsimple = simple(Bfull);
    Bsols{i} = Bsimple;
    Bsubs = subs(Bsimple, knowns, known_vals);
    for r = 1:3
        for c = 1:4
            numvar = numel(symvar(Bsubs(r,c)));
            numvars{i}(r,c) = numvar;
        end
    end
    numvars{i}
end

end


if 0

Beqs = {};
% BA * B06 * BB = BC
for i = 0:6
    for j = i:6
        BC = sym(eye(4));
        for k = i+1:j
            BC = BC * Bmats{k};
        end
        BA = sym(eye(4));
        for k = i:-1:1
            BA = BA * Bmats{k}^-1;
        end
        BB = sym(eye(4));
        for k = 6:-1:j+1
            BB = BB * Bmats{k}^-1;
        end
        Bfull = BA * B06 * BB - BC
        Beqs{end+1} = Bfull;
    end
end

end

if 0

alleqs = sym([]);
for i = 1:6
    knownvars = i-1
    foundsol = 0;
    for j = 1:numel(Beqs)
        Bsimple = Beqs{j};
        %Bsubs = subs(Bsimple, knowns, known_vals);
        Bsubs = simplify(Bsimple, 'Seconds', 0.3)
        %Bsubs = simple(Bsubs);
        j
        for r = 1:3
            for c = 1:4
                symvars = symvar(Bsubs(r,c));
                unknowns = setdiff(symvars, knowns)
                numvar = numel(unknowns);
                if numvar == 1
                    sprintf('j: %d, r: %d, c:%d\n', j, r, c)
                    alleqs(end+1) = Beqs{j}(r,c)
                    symvars
                    knowns(end+1) = unknowns(1);
                    known_vals(end+1) = 0.5;
                    foundsol = 1;
                    break
                end
            end
            if foundsol
                break
            end
        end
        if foundsol
            break
        end
    end
    if ~foundsol
        nosolutionfound = 1
        break
    end
end
alleqs
alleqssimp = simplify(alleqs)

end

if 0

Bsimples = {};
eq_list = {};
symvars_list = {};
for j = 1:numel(Beqs)
    j
    Bsimple = simplify(Beqs{j}, 'Seconds', 0.3);
    Bsimples{j} = Bsimple;
    for r = 1:3
        for c = 1:4
            symvars = symvar(Bsimple(r,c));
            eq_list{end+1} = Bsimple(r,c);
            symvars_list{end+1} = symvars;
        end
    end
end

end

if 1

cur_knowns = [x1, x2, x3, y1, y2, y3, z1, z2, z3, p1, p2, p3, d1, a2, a3, d4, d5, d6];
elim_order = [q1, q5, q6, q3, q2, q4];
possib_eqs = {};
skip = 0;
for i = 1:skip
    cur_knowns(end+1) = elim_order(i);
end
for i = skip+1:6
    knownvars = i-1
    possib_eqs{i} = {};
    for j = 1:numel(symvars_list)
        symvars = symvars_list{j};
        unknowns = setdiff(symvars, cur_knowns)
        if numel(unknowns) == 1
            if unknowns(1) == elim_order(i)
                possib_eqs{i}{end+1} = eq_list{j};
                eq_list{j}
            end
        end
    end
    cur_knowns(end+1) = elim_order(i);
end

end

