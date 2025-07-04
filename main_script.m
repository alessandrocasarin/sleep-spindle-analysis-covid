%% Progetto 1 - Neurobotics & Neurorehabilitation
%
%

clc 
close all
clearvars

%% DATA LOADING 

fs = 250;                       % Freqenza di campionnamento
load("Data/CTRL033_nap.mat");
ctrl_eeg = EEG;                 % Dati sogetto CONTROLLO
clear("EEG")
load("Data/ICU023_nap.mat");
icu_eeg = EEG;                  % Dati soggetto in ICU
clear("EEG")

%Vettori Tempo 
t_ctrl = (0 : 1/fs : (size(ctrl_eeg,1)-1)/fs);
t_icu = (0 : 1/fs : (size(icu_eeg,1)-1)/fs);

%% DATA FILERING for CONTROL
% Filtraggio del segnale controllo
% nelle bande lente (9-12)Hz e veloci (12-16)Hz

% Filtraggio causale con ordine 4

% Slow Tracks 9-12 Hz

ctrl_slow = []; %Inizializzo vettore

fc = 9;
[b,a] = butter(4,fc/(fs/2),"high"); % butterworth filter
ctrl_slow = filtfilt(b,a,ctrl_eeg); % forward-backward filtering 

fc = 12;
[b,a] = butter(4,fc/(fs/2),"low");
ctrl_slow = filtfilt(b,a, ctrl_slow);

clear("a","b")

% Fast Tracks 12-16 Hz

ctrl_fast = [];

fc = 12;
[b,a] = butter(4,fc/(fs/2),"high");
ctrl_fast = filtfilt(b,a,ctrl_eeg);

fc = 16;
[b,a] = butter(4,fc/(fs/2),"low");
ctrl_fast = filtfilt(b,a, ctrl_fast);

clear("a","b", "fc")

%% DATA FILERING for ICU

% Filtraggio del segnale icu
% nelle bande lente (9-12)Hz e veloci (12-16)Hz

% Filtraggio causale con ordine 4


% Slow Tracks 9-12 Hz

icu_slow = []; % Inizializzo il vettore

fc = 9;
[b,a] = butter(4,fc/(fs/2),"high"); % Butterworth filtering 
icu_slow = filtfilt(b,a,icu_eeg); % Forward-backward filtering

fc = 12;
[b,a] = butter(4,fc/(fs/2),"low");
icu_slow = filtfilt(b,a, icu_slow);

clear("a","b")

% Fast Tracks 12-16 Hz

icu_fast = [];

fc = 12;
[b,a] = butter(4,fc/(fs/2),"high");
icu_fast = filtfilt(b,a,icu_eeg);

fc = 16;
[b,a] = butter(4,fc/(fs/2),"low");
icu_fast = filtfilt(b,a, icu_fast);

clear("a","b", "fc")

%% Spindles Identification for CONTROL

spindle_wnd = 0.5*fs; % Ampiezza dello spindle (durata*fs)

% Load dei dati di posizione degli spindle
load("Data/spindles_timing_033.mat");
ctrl_fast_timing = fast;
ctrl_slow_timing = slow;
clear("fast","slow")

% Slow Spindles

ctrl_slow_spindles = {}; % Inizializzazione struct

%Selezionare le porzioni di segnale relative agli spindle
for spindle=1:length(ctrl_slow_timing)
    start_sample = round(ctrl_slow_timing(spindle,1)); % Campione iniziale
    end_sample = start_sample + spindle_wnd - 1;            % Campione finale
    ctrl_slow_spindles{end + 1} = ctrl_slow(start_sample:end_sample,:);
end
CTRL_SLOW_SPINDLES = cat(3, ctrl_slow_spindles{:}); % Concatenzazione 3D
clear("ctrl_slow_spindles", "start_sample", "end_sample", "spindle")

% Fast Spindles

ctrl_fast_spindles = {};

for spindle=1:length(ctrl_fast_timing)
    start_sample = round(ctrl_fast_timing(spindle,1));
    end_sample = start_sample + spindle_wnd - 1;
    ctrl_fast_spindles{end + 1} = ctrl_fast(start_sample:end_sample,:);

end
CTRL_FAST_SPINDLES = cat(3, ctrl_fast_spindles{:});
clear("ctrl_fast_spindles", "start_sample", "end_sample", "spindle")

%% Spindles Identification for ICU

spindle_wnd = 0.5*fs;

load("Data/spindles_timing_023.mat");
icu_fast_timing = fast;
icu_slow_timing = slow;
clear("fast","slow")

% Slow Spindles

icu_slow_spindles = {};

for spindle=1:length(icu_slow_timing)
    start_sample = round(icu_slow_timing(spindle,1));
    end_sample = start_sample + spindle_wnd - 1;
    icu_slow_spindles{end + 1} = icu_slow(start_sample:end_sample,:);

end
ICU_SLOW_SPINDLES = cat(3, icu_slow_spindles{:});
clear("icu_slow_spindles", "start_sample", "end_sample", "spindle")

% Fast Spindles

icu_fast_spindles = {};

for spindle=1:length(icu_fast_timing)
    start_sample = round(icu_fast_timing(spindle,1));
    end_sample = start_sample + spindle_wnd - 1;
    icu_fast_spindles{end + 1} = icu_fast(start_sample:end_sample,:);

end
ICU_FAST_SPINDLES = cat(3, icu_fast_spindles{:});
clear("icu_fast_spindles", "start_sample", "end_sample", "spindle")

%% Channel Averaging CONTROL
% Media degli spindle(trial) secondo i canali
% Dimensione Matrice risultante : Ns * 1 * N_trial

% Slow Spindles

CTRL_SLOW_SPINDLES_AVG = zeros(size(CTRL_SLOW_SPINDLES,1), ...
1,size(CTRL_SLOW_SPINDLES,3)); % Inizializzo la matrice

spindle_count = size(CTRL_SLOW_SPINDLES,3); % Numero di spindle

for spindle=1:spindle_count
    data_selection = CTRL_SLOW_SPINDLES(:,:,spindle); % Singolo trial
    channels = size(CTRL_SLOW_SPINDLES,2); % Numero canali del trial
    result = zeros(size(CTRL_SLOW_SPINDLES,1),1); % Risultato temporaneo
    for chn=1:channels
        result = result + data_selection(:,chn);
    end
    result = result ./ channels; % Calcolo della media lungo i canali
    CTRL_SLOW_SPINDLES_AVG(:,1,spindle) = result;
end

clear("spindle", "spindle_count", "data_selection", "channels", "result", "chn")

% Fast Spindles

CTRL_FAST_SPINDLES_AVG = zeros(size(CTRL_FAST_SPINDLES,1), ...
    1,size(CTRL_FAST_SPINDLES,3));

spindle_count = size(CTRL_FAST_SPINDLES,3);

for spindle=1:spindle_count
    data_selection = CTRL_FAST_SPINDLES(:,:,spindle);
    channels = size(CTRL_FAST_SPINDLES,2);
    result = zeros(size(CTRL_FAST_SPINDLES,1),1);
    for chn=1:channels
        result = result + data_selection(:,chn);
    end
    result = result ./ channels;
    CTRL_FAST_SPINDLES_AVG(:,1,spindle) = result;
end
clear("spindle", "spindle_count", "data_selection", "channels", "result", "chn")

%% Channel Averaging ICU

% Slow Spindles

ICU_SLOW_SPINDLES_AVG = zeros(size(ICU_SLOW_SPINDLES,1), ...
    1,size(ICU_SLOW_SPINDLES,3));

spindle_count = size(ICU_SLOW_SPINDLES,3);

for spindle=1:spindle_count
    data_selection = ICU_SLOW_SPINDLES(:,:,spindle);
    channels = size(ICU_SLOW_SPINDLES,2);
    result = zeros(size(ICU_SLOW_SPINDLES,1),1);
    for chn=1:channels
        result = result + data_selection(:,chn);
    end
    result = result ./ channels;
    ICU_SLOW_SPINDLES_AVG(:,1,spindle) = result;
end
clear("spindle", "spindle_count", "data_selection", "channels", "result", "chn")

% Fast Spindles

ICU_FAST_SPINDLES_AVG = zeros(size(ICU_FAST_SPINDLES,1), ...
    1,size(ICU_FAST_SPINDLES,3));

spindle_count = size(ICU_FAST_SPINDLES,3);

for spindle=1:spindle_count
    data_selection = ICU_FAST_SPINDLES(:,:,spindle);
    channels = size(ICU_FAST_SPINDLES,2);
    result = zeros(size(ICU_FAST_SPINDLES,1),1);
    for chn=1:channels
        result = result + data_selection(:,chn);
    end
    result = result ./ channels;
    ICU_FAST_SPINDLES_AVG(:,1,spindle) = result;
end
clear("spindle", "spindle_count", "data_selection", "channels", "result", "chn")
clear("max_peak","peak_loc", "peaks_pos","peaks")

%% SPECTRA for CONTROL
% Calcolo degli spettri in frequenza tramite funzione p-welch degli
% spindles

% Slow Spindles

ctrl_slow_spindle_spectra = {}; % Inizializzazione Struct
ctrl_slow_spectra_select  = {}; % Inizializzazione struct degli spettri 
                                % ammissibili

spindle_count = size(CTRL_SLOW_SPINDLES_AVG,3); % Numero degli spindle

frequency = (0.6:0.1:20); % Frequenze di interesse 

figure()
hold on
for spindle=1:spindle_count % Per ogni spindle
    data_selection = CTRL_SLOW_SPINDLES_AVG(:,1,spindle); % Singolo trial

    % Calcolo dello spettro
    [pxx, f] = pwelch(data_selection,length(data_selection),[],frequency,fs);
   
    [peaks, peaks_pos] = findpeaks(pxx); % Individuo i punti di massimo
    max_peak = find(peaks == max(peaks)); % Seleziono il massimo dei massimi
    peak_loc = peaks_pos(max_peak)*(f(2)-f(1)); % Calcolo frequenza del massimo picco

    % Range di frequenze accettabili
    f_min = 8;
    f_max = 12;

    % Se il picco ricade nel range si salva lo spindle (blue) altrimenti si
    % scarta (rosso)
    if peak_loc >= f_min && peak_loc <= f_max
        plot(f,pxx, 'b-');
        ctrl_slow_spectra_select{end + 1} = CTRL_SLOW_SPINDLES(:,:,spindle);
    else
        plot(f, pxx, 'r-')
    end
    grid on
    ctrl_slow_spindle_spectra{end + 1} = pxx;
end
xline(9,'r--') % Visualizza limite frequenze
xline(12, 'r--')
xlabel('freq. Hz')
ylabel('Amplitude')
title('CONTROL: Slow Spindle SPECTRAS')
hold off

%Calcolo delle matrici 3D
CTRL_SLOW_SPINDLE_SPECTRA = cat(3,ctrl_slow_spindle_spectra{:});
CTRL_SLOW_SPINDLE_SELECT = cat(3, ctrl_slow_spectra_select{:});
clear("ctrl_slow_spindle_spectra", "ctrl_slow_spectra_select", "spindle_count", "spindle", "data_selection","pxx")
clear("max_peak","peak_loc", "peaks_pos","peaks")

% Fast Spindles

ctrl_fast_spindle_spectra = {};
ctrl_fast_spectra_select  = {};

spindle_count = size(CTRL_FAST_SPINDLES_AVG,3);

frequency = (0.6:0.1:20);

figure()
hold on
for spindle=1:spindle_count
    data_selection = CTRL_FAST_SPINDLES_AVG(:,1,spindle);
    [pxx, f] = pwelch(data_selection,length(data_selection),[],frequency,fs);
    [peaks, peaks_pos] = findpeaks(pxx);
    max_peak = find(peaks == max(peaks));
    peak_loc = peaks_pos(max_peak)*(f(2)-f(1));

    f_min = 12;
    f_max = 16;

    if peak_loc >= f_min && peak_loc <= f_max
        plot(f,pxx, 'b-');
        ctrl_fast_spectra_select{end + 1} = CTRL_FAST_SPINDLES(:,:,spindle);
    else
        plot(f, pxx, 'r-')
    end
    grid on
    ctrl_fast_spindle_spectra{end + 1} = pxx;
end
xline(12,'r--')
xline(16, 'r--')
xlabel('freq. Hz')
ylabel('Amplitude')
title('CONTROL: Fast Spindle SPECTRAS')
hold off
CTRL_FAST_SPINDLE_SPECTRA = cat(3,ctrl_fast_spindle_spectra{:});
CTRL_FAST_SPINDLE_SELECT = cat(3, ctrl_fast_spectra_select{:});
clear("ctrl_fast_spindle_spectra", "ctrl_fast_spectra_select","spindle_count", "spindle", "data_selection","pxx")
clear("max_peak","peak_loc", "peaks_pos","peaks")

%% SPECTRA for ICU

% Slow Spindles

icu_slow_spindle_spectra  = {};
icu_slow_spectra_select = {};

spindle_count = size(ICU_SLOW_SPINDLES_AVG,3);

frequency = (0.6:0.1:20);

figure()
hold on
for spindle=1:spindle_count
    data_selection = ICU_SLOW_SPINDLES_AVG(:,1,spindle);
    [pxx, f] = pwelch(data_selection,length(data_selection),[],frequency,fs);
    [peaks, peaks_pos] = findpeaks(pxx);
    max_peak = find(peaks == max(peaks));
    peak_loc = peaks_pos(max_peak)*(f(2)-f(1));

    f_min = 8;
    f_max = 12;

    if peak_loc >= f_min && peak_loc <= f_max
        plot(f,pxx, 'b-');
        icu_slow_spectra_select{end + 1} = ICU_SLOW_SPINDLES(:,:,spindle);
    else
        plot(f, pxx, 'r-')
    end
    grid on
    icu_slow_spindle_spectra{end + 1} = pxx;
end
xline(9,'r--')
xline(12, 'r--')
xlabel('freq. Hz')
ylabel('Amplitude')
title('ICU: Slow Spindle SPECTRAS')
hold off
ICU_SLOW_SPINDLE_SPECTRA = cat(3,icu_slow_spindle_spectra{:});
ICU_SLOW_SPINDLE_SELECT = cat(3,icu_slow_spectra_select{:});
clear("ctrl_slow_spindle_spectra","icu_slow_spectra_select", "spindle_count", "spindle", "data_selection","pxx")
clear("max_peak","peak_loc", "peaks_pos","peaks")

% Fast Spindles

icu_fast_spindle_spectra = {};
icu_fast_spectra_select = {};

spindle_count = size(ICU_FAST_SPINDLES_AVG,3);

frequency = (0.6:0.1:20);

figure()
hold on
for spindle=1:spindle_count
    data_selection = ICU_FAST_SPINDLES_AVG(:,1,spindle);
    [pxx, f] = pwelch(data_selection,length(data_selection),[],frequency,fs);
    [peaks, peaks_pos] = findpeaks(pxx);
    max_peak = find(peaks == max(peaks));
    peak_loc = peaks_pos(max_peak)*(f(2)-f(1));

    f_min = 12;
    f_max = 16;

    if peak_loc >= f_min && peak_loc <= f_max
        plot(f,pxx, 'b-');
        icu_fast_spectra_select{end + 1} = ICU_FAST_SPINDLES(:,:,spindle);
    else
        plot(f, pxx, 'r-')
    end
    grid on
    icu_fast_spindle_spectra{end + 1} = pxx;
end
xline(12,'r--')
xline(16, 'r--')
xlabel('freq. Hz')
ylabel('Amplitude')
title('ICU: Fast Spindle SPECTRAS')
hold off
ICU_FAST_SPINDLE_SPECTRA = cat(3,icu_fast_spindle_spectra{:});
ICU_FAST_SPINDLE_SELECT = cat(3,icu_fast_spectra_select{:});
clear("icu_fast_spindle_spectra","icu_fast_spectra_select","spindle_count", "spindle", "data_selection","pxx")
clear("max_peak","peak_loc", "peaks_pos","peaks")

%% SPECTRAL AVERAGING - CONTROL
% Media degli spindle lungo i trial selezionati
% Dimensione della matrice risultante: Ns*Nch

% Slow Spindles

%Inizializzo la matrice
CTRL_SLOW_SPINDLE_MEAN = zeros(size(CTRL_SLOW_SPINDLE_SELECT,1),size(CTRL_SLOW_SPINDLE_SELECT,2));

spindle_count = size(CTRL_SLOW_SPINDLE_SELECT,3); % Numero di trial selezionati

%Calcolo la media
for spindle=1:spindle_count
    CTRL_SLOW_SPINDLE_MEAN(:,:) = CTRL_SLOW_SPINDLE_MEAN + CTRL_SLOW_SPINDLE_SELECT(:,:,spindle);
end
CTRL_SLOW_SPINDLE_MEAN = CTRL_SLOW_SPINDLE_MEAN ./ spindle_count;

% Fast Spindles

CTRL_FAST_SPINDLE_MEAN = zeros(size(CTRL_FAST_SPINDLE_SELECT,1),size(CTRL_FAST_SPINDLE_SELECT,2));

spindle_count = size(CTRL_FAST_SPINDLE_SELECT,3);

for spindle=1:spindle_count
    CTRL_FAST_SPINDLE_MEAN(:,:) = CTRL_FAST_SPINDLE_MEAN + CTRL_FAST_SPINDLE_SELECT(:,:,spindle);
end
CTRL_FAST_SPINDLE_MEAN = CTRL_FAST_SPINDLE_MEAN ./ spindle_count;

clear("spindle", "spindle_count")

%% SPECTRAL AVERAGING - ICU

% Slow Spindles

ICU_SLOW_SPINDLE_MEAN = zeros(size(ICU_SLOW_SPINDLE_SELECT,1),size(ICU_SLOW_SPINDLE_SELECT,2));

spindle_count = size(ICU_SLOW_SPINDLE_SELECT,3);

for spindle=1:spindle_count
    ICU_SLOW_SPINDLE_MEAN(:,:) = ICU_SLOW_SPINDLE_MEAN + ICU_SLOW_SPINDLE_SELECT(:,:,spindle);
end
ICU_SLOW_SPINDLE_MEAN = ICU_SLOW_SPINDLE_MEAN ./ spindle_count;

% Fast Spindles

ICU_FAST_SPINDLE_MEAN = zeros(size(ICU_FAST_SPINDLE_SELECT,1),size(ICU_FAST_SPINDLE_SELECT,2));

spindle_count = size(ICU_FAST_SPINDLE_SELECT,3);

for spindle=1:spindle_count
    ICU_FAST_SPINDLE_MEAN(:,:) = ICU_FAST_SPINDLE_MEAN + ICU_FAST_SPINDLE_SELECT(:,:,spindle);
end
ICU_FAST_SPINDLE_MEAN = ICU_FAST_SPINDLE_MEAN ./ spindle_count;

clear("spindle", "spindle_count")

%% TIME AVERAGING
% Media dei risultati lungo il tempo
% Dimensione della matrice risultante: 1*Nch

% Control - Slow

CTRL_SLOW_EEG = zeros(1,size(CTRL_SLOW_SPINDLE_MEAN,2));

t_count = size(CTRL_SLOW_SPINDLE_MEAN,1);
for t=1:t_count
    CTRL_SLOW_EEG = CTRL_SLOW_EEG + abs(CTRL_SLOW_SPINDLE_MEAN(t,:));
end
CTRL_SLOW_EEG = CTRL_SLOW_EEG ./ t_count;

% Control - Fast

CTRL_FAST_EEG = zeros(1,size(CTRL_FAST_SPINDLE_MEAN,2));

t_count = size(CTRL_FAST_SPINDLE_MEAN,1);
for t=1:t_count
    CTRL_FAST_EEG = CTRL_FAST_EEG + abs(CTRL_FAST_SPINDLE_MEAN(t,:));
end
CTRL_FAST_EEG = CTRL_FAST_EEG ./ t_count;

% ICU - Slow

ICU_SLOW_EEG = zeros(1,size(ICU_SLOW_SPINDLE_MEAN,2));

t_count = size(ICU_SLOW_SPINDLE_MEAN,1);
for t=1:t_count
    ICU_SLOW_EEG = ICU_SLOW_EEG + abs(ICU_SLOW_SPINDLE_MEAN(t,:));
end
ICU_SLOW_EEG = ICU_SLOW_EEG ./ t_count;

% ICU - Fast

ICU_FAST_EEG = zeros(1,size(ICU_FAST_SPINDLE_MEAN,2));

t_count = size(ICU_FAST_SPINDLE_MEAN,1);
for t=1:t_count
    ICU_FAST_EEG = ICU_FAST_EEG + abs(ICU_FAST_SPINDLE_MEAN(t,:));
end
ICU_FAST_EEG = ICU_FAST_EEG ./ t_count;

clear("t", "t_count")

%% TOPOPLOT 
% Visualizzazione in eeglab delle topografie dei 4 casi:
% CONTROLLO: Slow & Fast
% ICU: Slow & Fast

% Control - Slow
addpath('C:\eeglab2022.1')
eeglab % Apro eeglab (da aggiungere a matlabpath)

figure()
EEG.chanlocs = readlocs('Data/GSN_204.sfp'); % Carico file channel locations
subplot(221)
topoplot(CTRL_SLOW_EEG',EEG.chanlocs, 'style', 'both'); % visualizzo topoplot
colormap("jet")
colorbar
caxis([0 max(CTRL_SLOW_EEG')]);
title("CTRL - Slow Spindles")
subplot(222)
topoplot(CTRL_FAST_EEG',EEG.chanlocs, 'style', 'both');
colormap("jet")
colorbar
caxis([0 max(CTRL_FAST_EEG')]);
title("CTRL - Fast Spindles")
subplot(223)
topoplot(ICU_SLOW_EEG',EEG.chanlocs, 'style', 'both');
colormap("jet")
colorbar
caxis([0 max(ICU_SLOW_EEG')]);
title("ICU - Slow Spindles")
subplot(224)
topoplot(ICU_FAST_EEG',EEG.chanlocs, 'style', 'both');
colormap("jet")
colorbar
caxis([0 max(ICU_FAST_EEG')]);
title("ICU - Fast Spindles")