function electroDerivada = funcionDerivar(electro, Fs)
    num_der = [1 2 0 -2 -1];
    num_der = (Fs/8) * num_der;
    den_der = 1;
    [Hder,wder] = freqz(num_der,den_der);
    semilogx(wder* (Fs / (2 * pi)),20*log10(abs(Hder)))
    title('Derivador')
    xlabel('Frecuencia')
    ylabel('Amplitud')
  
    electroDerivada = filter(num_der, den_der, electro);
end