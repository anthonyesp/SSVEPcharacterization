# SSVEP-BCI characterization

## short description
Data and code for the metrological characterization of a low-cost wearable brain-computer interface relying on SSVEP.

Data are provided as .mat files.
Scripts were implemented in Matlab R2020b, but they could be used with previous versions of Matlab. 

## less short description
The low-cost wearable brain-computer interface operation relies on steady-state visually evoked potentials (SSVEP).
This system was implemented with:
- augmented reality glasses Espon Moverio BT-200, used for flickering stimuli (square white icons) presentation;
- electroencephalograph Olimex EEG-SMT with dry electrodes.


These devices were characterized independently of the human user.
Hence the uploaded data are related to:
- measures of light intensity associated with the flickering stimuli;
- measures of voltage (microvolts) done with the low-cost electroencephalograph;


The provided code allows to:
- characterize the flickering stimuli in terms of spectral content (amplitude) - see 'BT200_characterization_main.m';
- characterize the low-cost electroencephalograph in terms of linearity and gain error - see 'EEG_SMT_characterization_main.m';


## how to use
In order to reproduce our results as a starting point for your analyses, you can simply follow these steps:
1. download data and scripts (code --> download zip)
2. extract zip folder on your PC
3. for an example of flickering stimuli characterization, open "BT200_characterization_main.m" in the "BT-200 characterization" folder, and be sure to have the "BT200 frequency scan" subfolder with data (.mat files)
4. for an example of electroencephalograph characterization, open "EEG_SMT_characterization_main.m" in the "EEG-SMT characterization" folder, and be sure to have the "frequency scan" and "linearity" subfolders with data (.mat files); the subfolder "sinefit" is provided as well for the sine fitting algorithm needed within the script
5. enjoy

## further details

1. P. Arpaia, L. Callegaro, A. Cultrera, A. Esposito, and M. Ortolano, "Metrological characterization of a low-cost electroencephalograph for wearable neural interfaces in industry 4.0 applications", 2021 IEEE International Workshop on Metrology for Industry 4.0 & IoT (MetroInd4.0&IoT), [doi: 10.1109/MetroInd4.0IoT51437.2021.9488445](https://ieeexplore.ieee.org/document/9488445)
2. A. Cultrera, P. Arpaia, L. Callegaro, A. Esposito, and M. Ortolano, "Smart glasses and visually evoked potentials applications: characterisation of the optical output for different display technologies", MDPI Engineering Proceedings, [doi: 10.3390/ecsa-8-11263](https://www.mdpi.com/2673-4591/10/1/33/htm)
3. P. Arpaia, L. Callegaro, A. Cultrera, A. Esposito, and M. Ortolano, "Metrological characterization of consumer-grade equipment for wearable brain-computer interfaces and extended reality", IEEE Transactions on Instrumentation and Measurement, [doi: 10.1109/TIM.2021.3127650](https://ieeexplore.ieee.org/document/9612173)
4. "Reactive brain-computer interfaces and augmented reality project" on [ResearchGate](https://www.researchgate.net/project/Reactive-brain-computer-interfaces-and-augmented-reality)
