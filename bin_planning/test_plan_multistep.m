
t_duration = 400; % number of seconds in distributions
rate = 10; % time points / sec

% gaussian distributions for bin steps
% the last parameter weights the distribution by the probability
% the person is on this branch
%           beg mean, beg sig, end mean, end sig, bin prob,
% (will only display 50% of t_duration in visualization)
distribs = [       5,       4,       15,      4,      1.00;
                  17,       4,       19,      4,      1.00;
                  21,       4,       40,      4,      1.00;
                  40,       4,       45,      4,      1.00;
                  50,       4,       70,      4,      1.00;
                  75,       4,       95,      4,      1.00;
                 100,       4,      120,      4,      1.00;
];
bins = 1:7; % bin IDs
slot_states = [1, 2, 3]; % state of workspace slots (0 if empty, >0 if bin ID occupies)
nowtimesec = 20; % time (s) of current time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
probs = gen_data(distribs, rate, t_duration);
action = multistep(probs, bins, slot_states, nowtimesec, rate);
