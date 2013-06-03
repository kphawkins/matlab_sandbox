function [cost] = late_cost(t, td, startprobs, endprobs, binprob, nowtime)

if td > numel(t)
    cost = 1e10;
else
    if td >= nowtime+1
        costnow = sum(endprobs(nowtime+1:end)) * sum((t(td)-t(nowtime)).^2 .* startprobs(1:nowtime));
        costlater = sum((t(td)-t(nowtime+1:td)).^2 .* startprobs(nowtime+1:td));
        cost = binprob*(costnow + costlater);
    else
        cost = 0;
    end
end
