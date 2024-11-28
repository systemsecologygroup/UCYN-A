
clear all
clc
close all

%% Load files

load L_range.mat % L_range
load Nfix_1.mat % Nfix_1
load Nfix.mat % Nfix

%% Calculate max and min for the range of N2 fixation rate

for i = 1:length(L_range)
    Nfix_1_min(i) = min(min(Nfix_1(i,:,:)));
    Nfix_1_max(i) = max(max(Nfix_1(i,:,:)));
end

Nfix_1_max = smooth(L_range,Nfix_1_max,0.15,'loess');
Nfix_1_max = Nfix_1_max';

%% Plot the range of N2 fixation rate 

patch([L_range fliplr(L_range)], [Nfix_1_min fliplr(Nfix_1_max)], 'b','EdgeColor','none')
alpha(0.1)
hold on;

%% Plot N2 fixation rate 

plot(L_range,Nfix,'Color','b','LineWidth',2);
hold on;

%% Set plot requirements

xlim([0 800])
xlabel('PAR (\mumol m^{-2} s^{-1})')
xticks(0:200:800)
xticklabels({'0','200','400','600','800'})

ylim([0 10])
ylabel('N_2 fixation rate (fmol N cell^{-1} d^{-1})')
yticks(0:2:10)
yticklabels({'0','2','4','6','8','10'})

set(gca,'TickLength',[0.015, 0.01])
set(gca,'FontSize',13)

box on
hold off;

