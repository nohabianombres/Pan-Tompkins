function signalPasaAlta = filtro_pasalta(electro, Fs)
    %fil_5hz_num = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    %fil_5hz_den = [1 1];
    %[H5hz,w5hz] = freqz(fil_5hz_num,fil_5hz_den);

    %plot(w5hz* (Fs / (2 * pi)),20*log10(abs(H5hz)))
    %title('Pasa alta 5Hz')
    %xlabel('Frecuencia')
    %ylabel('Amplitud')



    Hd=finalfiltro5hz;
    [fil_5hz_num,fil_5hz_den] = sos2tf(Hd.sosMatrix, Hd.Scalevalue);

    [H5hz,w5hz] = freqz(fil_5hz_num,fil_5hz_den);
    figure
    plot(w5hz* (Fs / (2 * pi)),20*log10(abs(H5hz)))
    title('Pasa alta 11Hz')
    xlabel('Frecuencia')
    ylabel('Amplitud')
    signalPasaAlta = filter(fil_5hz_num, fil_5hz_den, electro);
end