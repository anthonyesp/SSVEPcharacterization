close all
clear
clc

% Circuit parameters
SH = 2;
FS = 4;
Nbit = 10;
QU = FS/2^Nbit;
    
% nominal gain from datasheet
GA = 7509.7;

% Vcal calculation from Vs
KRVD = 9.999317e-5;

% uncertainties
sKRVD = 13*10^-6;
tVS = 110*10^-6;            % obtained for FS = 1 V, reading ~ 1 V, no temperature coefficient

uBlin = 0.0091;             % relative uncertainty obtained from shifts in gain from linearity measures (uniform distribution)

format long                 % many digits are needed because of the Vs/VCAL precision

%% Linearity
% Consider all linearity files
files_linearity = dir('linearity\*.mat');
nf = length(files_linearity);

% VCAL rms
VCALrmsLIN = zeros(nf,1);

% Parameters from sinefit and uncertainty
sineParamLIN = zeros(nf,4);
u_vEEGLIN = zeros(nf,4);

% measured and fitted signal per each test
eegTotLIN = zeros(nf,2560);
eegFitLIN = zeros(nf,2560);

% number of samples for each test (it must be max 2560)
nSamples = zeros(nf,1);

% analyzing data by test frequency
freqLIN = 20;               % Hz
for k = 1:nf
    % Load data
    load(strcat(files_linearity(k,1).folder,'\',files_linearity(k,1).name))
    
    % scaling EEG codes
    eeg = 10^6*(eeg_CH1*QU-SH)/GA;

    % time and frequency
    nSaLIN = length(eeg);
    tLIN = 0:1/fs:(nSaLIN-1)/fs;  
    
    % sine fitting for more accurate measurement
    tw = nSaLIN/fs;
    w0 = 2*pi*freqLIN;
    tol = 2*pi*4/tw;
    
    % 4-parameters fitting
    [param4, eegFitLIN(k,1:nSaLIN), yerr] = sinefit(eeg,tLIN,[w0-tol, w0+tol]);
    
    % amplitude, frequency, offset, phase (rad)
    A = param4(2);
    O = param4(1);
    phi = param4(3);
    wu = param4(4);
    
    % fitting parameters
    % in order: [A O phi f]
    sineParamLIN(k,:) = [A, O, phi, round(wu/(2*pi),2)];
    
    % fitting parameters uncertainty
    sig = std(yerr);
    N = length(yerr);
    
    % in order: [u_A, u_O, u_phi, u_f]
    u_vEEGLIN(k,:) = [2*sig/sqrt(N) sig/sqrt(N) NaN NaN];
    
    % saving other useful data
    eegTotLIN(k,1:nSaLIN) = eeg;
    nSamples(k) = nSaLIN;
    VCALrmsLIN(k) = 10^6*V_CAL;
    
    clc
end

% rms values of the measured sinusoids
VEEGrmsLIN = sqrt(sineParamLIN(:,1).^2/2 + sineParamLIN(:,2).^2);

% associated uncertainties
u_VEEGrms = sqrt(((sineParamLIN(:,1).^2/2).*(u_vEEGLIN(:,1).^2) + (sineParamLIN(:,2).^2).*(u_vEEGLIN(:,2).^2))./(VEEGrmsLIN.^2));
kcov = 4;

% plotting linearity with an example of sinusoid as inset
figure
hold off
plot(VCALrmsLIN, VEEGrmsLIN,'k.-','MarkerFaceColor', 'k');
hold on
plot([0 140],[0 140],'k--');
errorbar(VCALrmsLIN, VEEGrmsLIN,kcov*u_VEEGrms,'k')
legend('experimental data','nominal response','Location','NorthWest')
xlabel('$V_\mathrm{CAL,rms}$ / $\mu$V', 'Interpreter','latex', 'FontSize', 12)
ylabel('$V_\mathrm{EEG,rms}$ / $\mu$V', 'Interpreter','latex', 'FontSize', 12)
set(gca, 'FontSize', 12)
xlim([0 120])
ylim([0 120])
ytickformat('%.0f')
grid
grid minor

k = 1;
t = 0:1/fs:(nSamples(k)-1)/fs;
axes('Position',[.6 .2 .25 .30])
box on
plot(t,eegTotLIN(k,1:nSamples(k)),'k')
xlim([0 0.2])
xticks([0 0.1 0.2])
ylim([-200 200])
grid
xtickformat('%.1f');
ytickformat('%.0f');
xlabel('t / s', 'Interpreter','latex', 'FontSize', 8)
ylabel('$V_\mathrm{EEG}$ / $\mu$V', 'Interpreter','latex', 'FontSize', 8)
set(gca, 'FontSize', 8)

% linear regression and test of linearity with ANOVA
N = length(VEEGrmsLIN);
b1 = (sum(VEEGrmsLIN.*VCALrmsLIN) - sum(VEEGrmsLIN)*sum(VCALrmsLIN)/N)/(sum(VCALrmsLIN.^2)-sum(VCALrmsLIN)^2/N);
yb = mean(VEEGrmsLIN);
xb = mean(VCALrmsLIN);
b0 = yb-b1*xb;

SST = sum((VEEGrmsLIN - yb).^2);
SSR = sum((b1*VCALrmsLIN+b0 - yb).^2);
SSE = sum((VEEGrmsLIN - (b1*VCALrmsLIN+b0)).^2);

gT = N - 1;
gE = N - 2;
gR = 1;

VART = SST/gT;
VARR = SSR/gR;
VARE = SSE/gE;

% Fisher's test for goodness of fitting
F = VARR/VARE;                              % figure of merit
Fb = finv(1-.05,gT,gE);                     % threshold F
p_value = 1-fcdf(F,gT,gE);

% r_lin = corr(VCALrmsLIN,VEEGrmsLIN);

%% Magnitude error (and frequency error)
% Consider all sine data files
files = dir('frequency scan\*.mat');
nf = length(files);

% test frequencies
freqs = zeros(nf,1);

% load one of the files to derive common acquisition parameters
load(strcat(files(1,1).folder,'\',files(1,1).name))
[~, MSGID] = lastwarn();

% number of samples (for EEG and Vs)
nSa_ees = length(eeg_CH1);

% measured Olimex EEG-SMT sampling frequency and time window
fs_ees = fs;
tw_ees = nSa_ees/fs_ees;

% time scales for the eess and REFs
t_ees = 0:1/fs_ees:(nSa_ees-1)/fs_ees;

% assuming the samples for Vs are fixed
nSa_sin = length(Vs);

% inizialize useful arrays
eesTot = zeros(nf,nSa_ees);                 % all acquired eegs
eesFit = zeros(nf,nSa_ees);                 % all fitted eegs
sinTot = zeros(nf,nSa_sin);                 % all acquired voltages
sinFit = zeros(nf,nSa_sin);                 % all reference voltages

% fitted parameters (amplitude, offset, phase, frequency)
eesParam = zeros(nf,4);
sinParam = zeros(nf,4);

% uncertainty of fitted parameters
u_eesParam = zeros(nf,4);
u_sinParam = zeros(nf,4);

warning('off', MSGID);

% analyzing different acquisitions (different test frequencies)
for k = 1:nf
    % load data per each acquisition
    load(strcat(files(k,1).folder,'\',files(k,1).name))
    freqs(k) = f_cal;
    
    % scaling EEG codes to microvolts
    ees = 10^6*(eeg_CH1*QU-SH)/GA;
    
    % saving EEG sinusoidal output for further analyses
    if (length(ees) ~= nSa_ees)
        error('the number of EEG samples is wrong');
    else
        eesTot(k,:) = ees;
    end

    % sine fitting for eegs (four-parameters)
    w0 = 2*pi*freqs(k);                       % nominal angular frequency
    tol = 0.02;                               % relative tolerance
    
    % amplitude, offset, phase (rad), frequency
    [param, eesFit(k,:), yerr] = sinefit(ees,t_ees,w0*[1-tol, 1+tol]);
    if (rms(yerr)/rms(ees) > 0.05)
        % repeat the fitting when the amplitude error is too high
        % start with the last fit frequency
        [param, eesFit(k,:), yerr] = sinefit(ees,t_ees,param(4)*[1-tol, 1+tol]);
        
        if (rms(yerr)/rms(ees) > 0.05)
            error(strcat('Fit error is too high even with 2 iterations for frequency ', num2str(f_cal)))
        end
    end
    
    A = param(2);
    O = param(1);
    phi = param(3);
    wu = param(4);
    
    % fitting parameters    
    eesParam(k,:) = [A, O, phi, wu/(2*pi)];
    
    % fitting parameters uncertainty
    sig = std(yerr);
    N = length(yerr);
    
    % in order: [u_A, u_O, u_phi, u_f]
    u_eesParam(k,:) = [2*sig/sqrt(N) sig/sqrt(N) NaN NaN];
    
    
    % time scale for Vs
    fs_sin = f_samp;
    tw_sin = nSa_sin/fs_sin;
    t_sin = 0:1/fs_sin:(nSa_sin-1)/fs_sin;
    
    % saving voltages for further analyses
    if (length(Vs) > nSa_sin)
        error('the number of VS samples is higher than expected');
    elseif (length(Vs) < nSa_sin)
        warning('the number of VS samples is lower than expected');
        sinTot(k,1:length(Vs)) = Vs;
        Vs = sinTot(k,:);           % adjust Vs to avoid error with sinefit
    else
        sinTot(k,:) = Vs;
    end
    
    % sine fitting for VSs (four-parameters)
    w0 = 2*pi*freqs(k);                     % nominal angular frequency is considered
    tol = 0.02;
    
    % amplitude, frequency, offset, phase (rad)
    [param, sinFit(k,:), yerr] = sinefit(Vs,t_sin,w0*[1-tol, 1+tol]);
    
    A = param(2); 
    O = param(1);
    phi = param(3);
    wu = param(4);
    
    % fitting parameters
    sinParam(k,:) = [A, O, phi, wu/(2*pi)];
    
    % fitting parameters uncertainty
    sig = std(yerr);
    N = length(yerr);
    
    % in order: [u_A, u_O, u_phi, u_f]
    u_sinParam(k,:) = [2*sig/sqrt(N) sig/sqrt(N) NaN NaN];
end

warning('on', MSGID);

% rms values with fitting
VEEGrmsAC = eesParam(:,1)/sqrt(2);
VEEGrmsDC = eesParam(:,2);
VEEGrms = sqrt(VEEGrmsAC.^2 + VEEGrmsDC.^2);

VSrmsAC = sinParam(:,1)/sqrt(2);
VSrmsDC = sinParam(:,2);
VSrms = sqrt(VSrmsAC.^2 + VSrmsDC.^2);

% evaluation of VCAL and its uncertainty
VCALrmsAC = 10^6*VSrmsAC*KRVD;
VCALrmsDC = 10^6*VSrmsDC*KRVD;
VCALrms = 10^6*VSrms*KRVD;


% uncertainties
u_VEEG_ac = u_eesParam(:,1);
u_VEEG_dc = u_eesParam(:,2);
u_VEEGrms = sqrt(((VEEGrmsAC.^2/2).*(u_VEEG_ac.^2) + (VEEGrmsDC.^2).*(u_VEEG_dc.^2))./(VEEGrms.^2));

u_VS_ac = u_sinParam(:,1);
u_VS_dc = u_sinParam(:,2);
u_VSrms = sqrt(((VSrmsAC.^2/2).*(u_VS_ac.^2) + (VSrmsDC.^2).*(u_VS_dc.^2))./(VSrms.^2));

% composite uncertainty VS (type A is sinefit, type B is HP3458A)
u_VSrmsC = sqrt(u_VSrms.^2 + tVS^2/3);

% u_VCAL_ac = VCALrmsAC.*sqrt(((u_VS_ac + tVS/sqrt(3))./VSrmsAC).^2 + sKRVD^2);           % recheck this formula
% u_VCAL_dc = VCALrmsDC.*sqrt(((u_VS_dc + tVS/sqrt(3))./VSrmsDC).^2 + sKRVD^2);           % recheck this formula
u_VCALrms = VCALrms.*sqrt((u_VSrmsC./VSrms).^2 + sKRVD^2);

% gain error plot and uncertainty
ERRrms = 100*(VEEGrms - VCALrms)./VCALrms;

% composite uncertainty VDUT (type A is sinefit, type B is from linearity gain shifts)
u_VEEGrmsC = sqrt(u_VEEGrms.^2 + (uBlin.*ERRrms).^2/3);

% coverage factor 4 --> 99.97% coverage for a Gaussian, in general at least 93.75%
kcov = 4;

uERRrms = 100*sqrt((u_VEEGrmsC.^2)./(VCALrms.^2) + ((VEEGrms./(VCALrms.^2)).^2).*(u_VCALrms.^2));

% Plotting results of tests with sinusoids
% clears the warning in the command window related to unattached EEG
clc

% Total rms error
figure
hold off
% hold on
plot(freqs,ERRrms, 'k-')
hold on
errorbar(freqs,ERRrms,kcov*uERRrms,'k.')
xlim([min(freqs)-1 max(freqs)+1])
% ylim([26 34])
xlabel('f / Hz', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\epsilon$ / \%', 'Interpreter','latex', 'FontSize', 12)
set(gca, 'FontSize', 12)
grid
grid minor

% frequency error
Err_f = (eesParam(:,4) - sinParam(:,4))./sinParam(:,4);
figure
hold off
% hold on
plot(freqs,100*Err_f, 'k-*')
xlim([min(freqs)-1 max(freqs)+1])
% ylim([26 34])
xlabel('$f_{CAL}$ / Hz', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\epsilon_f$ / \%', 'Interpreter','latex', 'FontSize', 12)
ytickformat('%.1f')
ylim([-0.5 0])
set(gca, 'FontSize', 12)
grid
grid minor