function signalPasaBaja = filtro_pasabaja (electro,Fs)
    fil_11hz_num = conv([1 0 0 0 0 0 -1], [1 0 0 0 0 0 -1]);
    fil_11hz_den = conv([1 -1], [1 -1]);
    [H11hz,w11hz] = freqz(fil_11hz_num,fil_11hz_den);
    figure
    plot(w11hz* (Fs / (2 * pi)),20*log10(abs(H11hz)))
    title('Pasa baja 5Hz')
    xlabel('Frecuencia')
    ylabel('Amplitud')
    signalPasaBaja = filter(fil_11hz_num, fil_11hz_den, electro);
end