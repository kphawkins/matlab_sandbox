%plan = create_plan(slot_states, deliv_seq)
%action_durs = ones(1,5) * 30;
%lower_bounds = round([0.1, action_durs] * rate);
%upper_bounds = ones(1,6) * inf;
%upper_bounds = [];
%guess = lower_bounds;
%guess = [10, 2, 2, 2, 2, 2] * rate;
%guess = lower_bounds * rate;
%opt_cost_fun(guess, slot_states, plan, t, probs, undodur, nowtime, 0)
%x_sol = fmincon(@(x) opt_cost_fun(x, slot_states, plan, t, probs, undodur, nowtime, 0), ...
%                guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
%best_times = cumsum(x_sol / rate)
%opt_cost_fun(x_sol, slot_states, plan, t, probs, undodur, nowtime, 0)
%opt_cost_fun([10, 5, 2, 2, 2, 2] * rate, slot_states, plan, t, probs, undodur, nowtime)
planning_params
deliv_seqs = generate_plans(t, beam_counts, probs, bins, slot_states, nowtime, endedweight, notbranchweight)

all_best_times = [];
all_costs = [];
all_plans = [];
for i = 1:size(deliv_seqs,1)
    plan = create_plan(slot_states, deliv_seqs(i,:));
    durations = create_durations(plan, traj_dur);
    lower_bounds = durations;
    lower_bounds(1) = lower_bounds(1) + nowtime;
    x_sol = fmincon(@(x) opt_cost_fun(x, slot_states, plan, t, probs, undodur, nowtime, 0), ...
                    lower_bounds, [], [], [], [], lower_bounds, [], [], options);
    best_times = cumsum(x_sol / rate);
    [cost, plan] = opt_cost_fun(x_sol, slot_states, plan, t, probs, undodur, nowtime, 1);
    deliver_sequence = deliv_seqs(i,:);
    all_best_times(i,:) = best_times;
    all_costs(i) = cost;
    all_plans(i,:) = plan;
end

[costs_sorted, cost_inds] = sort(all_costs);
for i = 1:size(deliv_seqs,1)
    i
    ind = cost_inds(i);
    cost = all_costs(ind)
    best_times = all_best_times(ind,:)
    durations = create_durations(plan, traj_dur);
    plan = all_plans(ind,:)
    figure(100+i)
    clf
    subplot(2,1,2)
    visualize_bin_probs(t, bins, probs, nowtime/rate, 100);
    subplot(2,1,1)
    title(sprintf('Cost: %.1f', cost))
    visualize_bin_activity(plan, [(best_times-durations/rate)', best_times'], bins, nowtime/rate, 100);
end
