%% Implementação e analise do filtro

num = zeros(1,33);
num(1) = -1; num(17) = 32; num(18)= -32 ; num(33) = 1;
den = 32*[1 -1];

% Resposta em frequência do filtro
freqz(num,den);

% Busca pela frequência de corte
[h,w_Norm] = freqz(num,den);
mag = 20*log10(abs(h));
w_piNorm = w_Norm/pi; % normalizada por pi
for k1 = 1:size(mag,1)
    if (mag(k1) >= -3) % Hp: '>=' & Lp: '<=' 
        break
    end
end

% k1 passa a ser o índice da primeira frequência menor que -3dB
cutoff_freq_w_norm = w_piNorm(k1);
cutoff_freq_norm = cutoff_freq_w_norm*pi/(2*pi); % Mudança de unidade de w para f com a relação w = 2*pi*f -> f = w/(2*pi)
cutoff_freq = cutoff_freq_norm*Fs/2;

plot(w_piNorm(k1),mag(k1),'*',w_piNorm,mag)
grid minor
ax = gca;
ax.YLim = [-100 20];
ax.XTick = 0:.5:2;
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')