    Fs = 200;
    num_der = [25 50 0 -50 -25];
    den_der = 1;
    [Hder,wder] = freqz(num_der,den_der);
    figure
    semilogx(wder* (Fs / (2 * pi)),20*log10(abs(Hder)))
    title('Derivador')
    xlabel('Frecuencia')
    ylabel('Amplitud')