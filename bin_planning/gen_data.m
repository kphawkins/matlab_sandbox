
% http://en.wikipedia.org/wiki/Normal-gamma_distribution
Tend = 200;
rate = 10;
sigalpha = 2;
meanvar = 4;
mumeans = [5, 15;
           17, 19;
           21, 40;
           40, 45;
           50, 70;
           75, 95];

varmean = 0;

sigbeta = sigalpha/meanvar;
N = Tend*rate;
numbins = size(mumeans,1);

t = linspace(0,Tend,N);

probs = cell(numbins,2);
for i = 1:numbins
    probs{i,1} = draw_normal(t, mumeans(i,1), varmean, sigalpha, sigbeta);
    probs{i,2} = draw_normal(t, mumeans(i,2), varmean, sigalpha, sigbeta);
end

figure(1)
clf
hold on
for i = 1:numbins
    plot(t,probs{i,1},'r')
    plot(t,probs{i,2},'b')
end
