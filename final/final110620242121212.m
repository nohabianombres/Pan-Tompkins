%% Parte a - Señal original

% Datos de entrada (señal de ejemplo)
%% Limpieza del entorno de trabajo
clear all;
clc;

%% Definir la ruta de la carpeta que contiene los archivos
folder_path = 'C:\Users\Jose\OneDrive - UCO\Desktop\SEMESTRE 6\PDS\final\archivos';  % Actualiza con la ruta correcta

%% Especificar el nombre del archivo MAT a cargar
mat_filename = '223m.mat';  % Actualiza con el nombre del archivo que deseas cargar

%% Construir el nombre completo del archivo MAT
mat_fullpath = fullfile(folder_path, mat_filename);

%% Verificar si el archivo MAT existe y cargar los datos
if exist(mat_fullpath, 'file') == 2
    % Cargar los datos del archivo MAT
    data_struct = load(mat_fullpath);
    
    % Suponiendo que la señal de ECG está almacenada en una variable llamada 'val'
    electro = data_struct.val; 
    
    % Mostrar los primeros 10 elementos de los datos
    fprintf('Primeros 10 elementos del archivo %s:\n', mat_filename);
    disp(electro(1:10));
    
    % Frecuencia de muestreo
    Fs = 360; % Actualiza esto con la frecuencia de muestreo correcta si es diferente
    
    L = length(electro);     % Longitud de la señal (numero de muestras)
    T = L/Fs;               % Duración de la señal
    t = linspace(0,T,L);
    
    figure;
    subplot(7, 1, 1);
    plot(electro);
    title('a) Señal Original');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Parte b - Salida del filtro pasa banda - cascada
    B_lp_1 = conv([1 0 0 0 0 0 -1], [1 0 0 0 0 0 -1]);
    A_lp_1 = [1 -2 1];
    electro_pb = filter(B_lp_1,A_lp_1,electro);
    electro_pb = electro_pb/max(electro_pb);

    %Highpass filter
    B_hp_1 = -[1 zeros(1,15) -32 zeros(1,15) -1];
    A_hp_1 = [1 1];
    electro_pb = filter(B_hp_1,A_hp_1,electro_pb);
    electro_pb = electro_pb/max(electro_pb);

    %Derivative filter
    %B_df_1 = [1 2 0 -2 -1]*Fs/8;
    %A_df_1 = [1];
    %ecg_df = filter(B_df_1,A_df_1,ecg_hp);
    %ecg_df  = ecg_df/max(ecg_df);
    % Filtro Pasa Bajas: Frecuencia de corte 12 Hz
    %fil_12hz_num = conv([1 0 0 0 0 0 -1], [1 0 0 0 0 0 -1]);
    %fil_12hz_den = conv([1 -1], [1 -1]);
    %electro_pb = filter(fil_12hz_num, fil_12hz_den, electro);
    
    % Filtro Pasa Altas: Frecuencia de corte 5 Hz
    %Hd = finalfiltro5hz();
    %[fil_5hz_num, fil_5hz_den] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);
    %electro_pb = filter(fil_5hz_num, fil_5hz_den, electro_pb);
    
    subplot(7, 1, 2);
    plot(electro_pb);
    title('b) Salida del Filtro Pasa Banda');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Parte c - Salida del diferenciador - Verificación de cambios en la pendiente de la señal del ECG
    
    % Filtro Derivador
    num_der = [45 90 0 -90 -45];
    den_der = 1;
    electro_der = filter(num_der, den_der, electro_pb);
    electro_der = electro_der/max(electro_der);
    
    subplot(7, 1, 3);
    plot(electro_der);
    title('c) Salida del Diferenciador');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Parte d - Cuadratura (se eleva al cuadrado bit a bit para resaltar los picos y reducir las Amp negativas)
    
    electro_sq = electro_der.^2;
    
    subplot(7, 1, 4);
    plot(electro_sq);
    title('d) Resultado del Proceso de Cuadratura');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Parte e - Integración Ventana móvil - Se suaviza la señal para identificar los complejos QRS
    
    N = 54; % Ventana de 150 ms (30 muestras a 200 Hz)
    num_samples = length(electro_sq);
    y = zeros(1, num_samples);
    for n = N:num_samples
        y(n) = (1 / N) * sum(electro_sq(n - (N - 1):n));
    end
    y(1:N-1) = 0;
    electro_mv = y;
    
    subplot(7, 1, 5);
    plot(electro_mv);
    title('e) Resultados de la Integración con Ventana Móvil');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Parte f - Señal de ECG original (tiene un retraso por el tiempo de procesamiento)
    delay = round(N / 2);
    electro_delayed = [zeros(1, delay), electro(1:end-delay)];
    
    subplot(7, 1, 6);
    plot(electro_delayed);
    title('f) Señal de ECG Original Retrasada');
    xlabel('Muestras');
    ylabel('Amplitud');
    
    %% Umbralización (Detección de QRS)
    
    thresholdSetF = 0.8 * max(electro_pb);  % Inicializando el umbral para la señal filtrada
    thresholdSetI = 0.8 * max(electro_mv);  % Inicializando el umbral para la señal integrada
    RRlow = 0;
    RRavg1 = 0;
    RRavg2 = 0;
    RR1 = [];
    SPKF = 0;
    NPKF = 0;
    SPKI = 0;
    NPKI = 0;
    refactoryPeriod = ceil(360/((Fs.^(-1))*1000));  % Periodo refractario de 200 ms
    prevQRSindex = 0;
    
    QRSs = [];
    

    % Encuentra los picos de los filtros usando findpeaks
    [~, peaksFiltidx] = findpeaks(electro_pb);
    
    % Inicializa la lista de valores con ceros
    peaksFiltro = zeros(size(electro_pb));
    
    % Asigna los valores de los picos del filtro
    peaksFiltro(peaksFiltidx) = electro_mv(peaksFiltidx);
    


    [~, peaksIntidx] = findpeaks(electro_pb);
    
    % Inicializa la lista de valores con ceros
    peaksIntegra = zeros(size(electro_mv));
    
    % Asigna los valores de los picos del filtro
    peaksIntegra(peaksIntidx) = electro_mv(peaksIntidx);

    RRmiss = 0;

    

    for i = 1:length(electro_mv) - 1
        % Detecta picos QRS en la señal filtrada
        if peaksFiltro(i) > thresholdSetF
            SPKF = (0.125 * peaksFiltro(i)) + 0.875 * SPKF;
            isQRSFilt = true; 
        else
            NPKF = (0.125 * peaksFiltro(i)) + 0.875 * NPKF;
            isQRSFilt = false;
        end
        
        % Detecta picos QRS en la señal integrada
        if peaksIntegra(i) > thresholdSetI
            SPKI = (0.125 * peaksIntegra(i)) + 0.875 * SPKI;
            isQRSInt = true;
        else
            NPKI = (0.125 * peaksIntegra(i)) + 0.875 * NPKI;
            isQRSInt = false;
        end
        
        % Determina si ambos picos (filtrada e integrada) coinciden y cumplen con el periodo refractario
        if isQRSFilt && isQRSInt && (i - prevQRSindex) > refactoryPeriod && (i - prevQRSindex) > RRlow
            % Busca el máximo pico en la ventana de detección alrededor de 'i'
            start_idx = max(1, i - round(refactoryPeriod / 2));
            end_idx = min(length(electro_delayed), i + round(refactoryPeriod / 2));
            
            % Find the index of the maximum value within the specified range
            [~, max_rel_idx] = max(electro_delayed(start_idx:end_idx));
            
            % Calculate the absolute index of the maximum value
            QRS_idx = start_idx + max_rel_idx - 1;
            
            % Append the detected QRS index to the QRSs list
            QRSs = [QRSs, QRS_idx];
            
            % Calcula el intervalo RR
            if prevQRSindex > 0
                RR = QRS_idx - prevQRSindex;
                RR1 = [RR1, RR];
                
                % Actualiza los promedios RR
                if length(RR1) >= 8
                    RRavg1 = mean(RR1(end-7:end));
                else
                    RRavg1 = mean(RR1);
                end
                
                RRavg2 = mean(RR1);
                
                % Umbrales de RR
                RRlow = 0.92 * RRavg2;
                RRhigh = 1.16 * RRavg2;
                RRmiss = 1.66 * RRavg2;
                
                % Ajusta los umbrales si se encuentra un intervalo RR grande (latido perdido)
                if RR > RRmiss
                    thresholdSetF = 0.5 * thresholdSetF;
                    thresholdSetI = 0.5 * thresholdSetI;
                end
            end
            
            prevQRSindex = QRS_idx;
        end
        
        % Actualiza los umbrales
        thresholdSetF = NPKF + 0.25 * (SPKF - NPKF);
        thresholdSetI = NPKI + 0.25 * (SPKI - NPKI);
    end
    
    subplot(7, 1, 7);
    plot(electro_delayed);
    hold on;
    plot(QRSs, electro_delayed(QRSs), 'ro');
    %title('g) Señal de ECG con Detección de QRS');
    xlabel('Muestras');
    ylabel('Amplitud');
    legend('ECG', 'QRS');



    %pulse_output = zeros(1, length(electro_delayed));
    %pulse_output(QRSs) = 1;
    
    %subplot(7, 1, 7);
    %stem(pulse_output);
    %title('g) Flujo de Pulsos de Salida');
    %xlabel('Muestras');
    %ylabel('Detección de QRS');

else
    fprintf('El archivo %s no se encuentra en la carpeta especificada.\n', mat_filename);
end
