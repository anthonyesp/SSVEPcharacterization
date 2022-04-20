close all
clear
clc

load('BT200_opt101_10Hz_test3.mat')
adc = 2.2*adc_CH1;

N = length(adc);
t = 0:1/fs:(N-1)/fs;

win = 'bh';

% windowing and coherent gain
if (strcmp(win,'flat'))
    WDW = transpose(flattopwin(N));
    CG = 0.22;
elseif (strcmp(win,'bh'))
    WDW = transpose(blackmanharris(N));
    CG = 0.42;
elseif (strcmp(win, 'hann'))
    WDW = transpose(hanning(N));
    CG = 0.5;
else
    % rect (no window)
    WDW = 1;
    CG = 1;
end

%     adc = movmean(adc,10);
ADC = abs(fft(WDW.*(adc-mean(adc))))/length(adc)/CG;
ADC2 = [ADC(1) 2*ADC(2:floor(end/2)) ADC(floor(end/2)+1)];

tw = N/fs;
f = 0:1/tw:fs-1/tw;
f2 = f(1:floor(end/2)+1);

% plot in time
figure
subplot(2,3,1)
plot(t,adc,'k')
xlim([0 1])
lim = 1.1*max(abs(adc));
%     ylim([-lim +lim])
ylim([0 3.0])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$t$ / s','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(a)')
grid

% plot in frequency
subplot(2,3,4)
semilogy(f2,ADC2,'k')
%     plot(f2,ADC2,'r')
xlim([0 100])
LIM = 1.1*max(ADC2);
%     ylim([0 LIM])
ylim([0.0001 1.0])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$f$ / Hz','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(d)')
grid
    

load('BT350_10Hz_1.mat')
adc = adc_CH1;
ADC = abs(fft(WDW.*(adc-mean(adc))))/length(adc)/CG;
ADC2 = [ADC(1) 2*ADC(2:floor(end/2)) ADC(floor(end/2)+1)];

% plot in time
subplot(2,3,2)
plot(t,adc,'k')
xlim([0 1])
lim = 1.1*max(abs(adc));
%     ylim([-lim +lim])
ylim([0 2.5])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$t$ / s','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(b)')
grid

% plot in frequency
subplot(2,3,5)
semilogy(f2,ADC2,'k')
%     plot(f2,ADC2,'r')
xlim([0 100])
LIM = 1.1*max(ADC2);
%     ylim([0 LIM])
ylim([0.0001 1.0])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$f$ / Hz','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(e)')
grid

load('holo_FLICK_10Hz_1.mat')
adc = adc_CH1;
ADC = abs(fft(WDW.*(adc-mean(adc))))/length(adc)/CG;
ADC2 = [ADC(1) 2*ADC(2:floor(end/2)) ADC(floor(end/2)+1)];

% plot in time
subplot(2,3,3)
plot(t,adc,'k')
xlim([0 1])
lim = 1.1*max(abs(adc));
%     ylim([-lim +lim])
ylim([0 3.5])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$t$ / s','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(c)')
grid

% plot in frequency
subplot(2,3,6)
semilogy(f2,ADC2,'k')
%     plot(f2,ADC2,'r')
xlim([0 100])
LIM = 1.1*max(ADC2);
%     ylim([0 LIM])
ylim([0.0001 1.0])
ytickformat('%.1f');
xtickformat('%.1f');
set(gca,'FontSize',18);
xlabel('$f$ / Hz','Interpreter','Latex','FontSize',18)
ylabel('$V_\mathrm{PT}$ / V','Interpreter','Latex','FontSize',18)
title('(f)')
grid