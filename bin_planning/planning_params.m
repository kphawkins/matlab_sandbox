slot_states = [1, 2, 3];
deliv_seq = 4:6;
%nowtime = 1;
nowtime = 20*rate;
beam_counts = [2, 2, 1];
traj_dur = 1*rate;
undodur = 2*traj_dur;
bins = 1:6;
endedweight = 10;
notbranchweight = 10;
options = optimset('Algorithm', 'active-set', 'FinDiffRelStep', 1, 'MaxFunEvals', 100);
