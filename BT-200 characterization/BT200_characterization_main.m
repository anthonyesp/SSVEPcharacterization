close all
clear
clc

% OPT101 parameters for voltage conversion
GOPT = 11e+6;                   % gain in V/A
ROPT = 0.45;                    % photodiode current vs input power @ 650 nm (A/W)
AOPT = 5.2e-2;                  % photodiode area in cm^2
uWcm2_to_lux = 1/(1.4641e-7);   % conversion factor from W/cm2 to lux

% conversion factor from ADC voltage to input lux
CF = uWcm2_to_lux/(GOPT*ROPT*AOPT);             

%% BT200 data analysis
% flickering data
load('BT200 frequency scan\BT200_opt101_10Hz_test1.mat')
BT200_time = adc_CH1;

% sampling frequency 'fs = 1000 Sa/s' is equal for all of them
% acquisition time should also be equal for all of them: tw = 10 s
N = length(BT200_time);
tw = N/fs;

% time scale
t = 0:1/fs:tw-1/fs;

% frequency scale
f = 0:1/tw:fs-1/tw;
f2 = f(1:floor(end/2)+1);

% windowing
WDWbh = transpose(blackmanharris(N));
CGbh = 0.42;

WDWha = transpose(hanning(N));
CGha = 0.50;

WDWno = 1;
CGno = 1.00;
    
% window to use
WDW = WDWha;
CG = CGha;

% spectra
TEMP = abs(fft(WDW.*(BT200_time-mean(BT200_time))))/N/CG;
BT200_freq = [TEMP(1) 2*TEMP(2:floor(end/2)) TEMP(floor(end/2)+1)];

[BT200pks, BT200locs] = findpeaks(BT200_freq,'MinPeakDistance',90);
BT200pks = BT200pks(2:16);
BT200locs = (BT200locs(2:16)-1)/tw;
if (BT200locs(end)/12 ~= 15)
    warning('unexpected last harmonic for BT200')
end
BT200_harm1 = BT200pks(1);

% SPECTRA COMPARISON (logarithmic plot)
figure
plot(t,BT200_time*CF,'k');
xlim([0 1])
ylim([0 35])
xlabel('t / s')
ylabel('luminous intensity / lux')
grid
xtickformat('%.1f')
ytickformat('%.0f')

figure
semilogy(f2,BT200_freq/BT200_harm1,'k')
xlim([0 100])
ylim([0 1.0])
hold on
xlabel('f / Hz')
ylabel('normalized harmonics')
grid
xtickformat('%.1f')
ytickformat('%.1f')

% harmonic ratio error log
harm_odd = [1 3 5 7 9 11 13 15];
harm_even = [1 0 3 0 5 0 7 0 9 0 11 0 13 0 15];

BT200ratios = BT200pks(1)./BT200pks;
BT200harm_err = 100*(BT200ratios - harm_even)./harm_even;

% linear frequency plot with harmonic ratios
figure
plot(f2,BT200_freq/BT200_harm1,'k')
xlim([0 100])
ylim([0 1.0])
hold on
xlabel('f / Hz')
ylabel('normalized harmonics')
grid
xtickformat('%.1f')
ytickformat('%.1f')
for i = 2:7
    text(BT200locs(i)+1.2,(1/BT200ratios(i))+0.04,strcat('1/',num2str(round(BT200ratios(i)),2)),'Color', 'k')
end
