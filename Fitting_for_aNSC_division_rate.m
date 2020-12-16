%% CURVE FITTING FOR gamma_a
%   This script calculates and plots a fitting curve for the decay in the fraction 
%   of aNSC tracks (experimentally observed. See figure S6 in Dray et al.) that did not divide by time t 



%% Table 1: Number of aNSC tracks dividing in a certain time range
Number_of_aNSC_tracks_that_divide_within_the_time_range = {'0 to 2 days';'3 to 5 days';'6 to 8 days';'9 to 11 days';'12 to 14 days';'16 days or more';'total'};
titi = [ 9 12 1 0 0 1 23]'; 
mimi = [ 6 10 7 3 2 1 29]';
bibi = [ 7 17 7 2 0 3 36]';
average = round(mean([titi mimi bibi],2),2,'significant'); 
STD = round(std([titi mimi bibi],0,2),2,'significant'); 
T1 = table(Number_of_aNSC_tracks_that_divide_within_the_time_range,titi,mimi,bibi,average,STD)

%% Table 2: fraction of aNSC tracks that did not divide by time t = (total number of aNSC tracks - Number of aNSC tracks that divided by time t)/total number of aNSC tracks
Number_of_aNSC_tracks_that_did_not_divide_by_time_t___ = {'t = 2 days';'t = 5 days';'t = 8 days';'t = 11 days';'t = 14 days';'t = 16 days'};
titi_remain(1,1) = titi(end,1) - titi(1,1);
mimi_remain(1,1) = mimi(end,1) - mimi(1,1);
bibi_remain(1,1) = bibi(end,1) - bibi(1,1);
for i = 2:length(titi)-1
    titi_remain(i,1) = titi_remain(i-1,1) - titi(i,1);
    mimi_remain(i,1) = mimi_remain(i-1,1) - mimi(i,1);
    bibi_remain(i,1) = bibi_remain(i-1,1) - bibi(i,1);
end
titi = titi_remain./titi(end,1);
mimi = mimi_remain./mimi(end,1);
bibi = bibi_remain./bibi(end,1);
average = round(mean([titi mimi bibi],2),2,'significant');
STD = round(std([titi mimi bibi],0,2),2,'significant');
T2 = table(Number_of_aNSC_tracks_that_did_not_divide_by_time_t___,titi,mimi,bibi,average,STD)

%% fitting
t = [ 0 2 5 8 11 14 16 ]';
f = [ 1 ; average];
f_stds = [ 0 ; STD ] ;
df = [10^4 ; 1./(STD(1:end-1).^2) ; 10^4];

Func = @(gamma_a,x)(exp(-gamma_a.*x));

[curve, goodness, output] = fit( t, f, Func, ...
    'StartPoint', 0.2, ...
    'Lower', 0, ...
    'Upper', 0.3, ...
    'Robust', 'Bisquare','Weights',df );
curve


%%  plot
figure('Renderer', 'painters', 'Position', [10 10 900 600])
p = plot( curve, t, f, 'bo');
set(p,'MarkerSize',8);
set(p,'linewidth',2);
hold on
e = errorbar(t,f,f_stds,'bo');
set(e,'linewidth',2);
set(e,'MarkerSize',0.01);
legend({'Data',' y=e^-^{\gamma_a}^t',' STD'},'FontSize',18)
xlabel('t [days]')
ylabel('fraction of aNSCs that did not divide')
title('Estimation of aNSCs division rate')
ax = gca;
ax.XTick = 0:1:16;
labels = string(ax.XTickLabel); % extract
labels(2:2:end) = nan; % remove every other one
ax.XTickLabel = labels; % set
ax.YTick = 0:0.2:1;
set(ax,'FontSize',18,'FontWeight','bold')