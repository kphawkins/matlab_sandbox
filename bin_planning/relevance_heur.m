function [h] = relevance_heur(t, startprobs, endprobs, binprob, nowtime, endedweight, notbranchweight)

probended = sum(endprobs(nowtime:end));
expectedstart = sum(startprobs .* t);
h = (t(nowtime)-expectedstart) + endedweight*(1-1/probended) + notbranchweight*(1-1/binprob);
