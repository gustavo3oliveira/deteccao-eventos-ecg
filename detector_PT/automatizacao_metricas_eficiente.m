% Author: Gustavo Henrique de Oliveira, Maria Ant�nia Marzo Barcel�, Roberto Costa Ceccato
% Date: 27/06/2021
% Description: Script principal para o projeto da Parte B (Prof. S�rgio)

% save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as vari�veis
% ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna cont�m sinal com Nt amostras j� calibradas em mV
% ts(Nt,1): [double] vetor com Nt instantes de tempo em s
% Fs: frequencia de amostragem em Hz
% ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
% anotacao ocorreu na amostra 270, entao p=ann(3) retornar� p=270.
% type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.      

close all 

%% Load do sinal no espa�o de trabalho
% save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as vari�veis
% ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna cont�m sinal com Nt amostras j� calibradas em mV
% ts(Nt,1): [double] vetor com Nt instantes de tempo em s
% Fs: frequencia de amostragem em Hz
% ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
%       anotacao ocorreu na amostra 270, entao p=ann(3) retornar� p=270.
% type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.  

n_sinais = [111:119]';
Sinal = n_sinais;

% Defini��o dos vetores que abrigar�o as m�tricas
%Sinal = zeros(size(n_sinais));
Media = zeros(size(n_sinais));
DP = zeros(size(n_sinais));
Media_norm = zeros(size(n_sinais));
DP_norm = zeros(size(n_sinais));
Nv  = zeros(size(n_sinais));
Nd = zeros(size(n_sinais));
FP_T = zeros(size(n_sinais));
p_FP = zeros(size(n_sinais));
FN_T = zeros(size(n_sinais));
p_FN = zeros(size(n_sinais));

% Relacionados ao IRR
Media_IRR = zeros(size(n_sinais));
DP_IRR = zeros(size(n_sinais));

for k0 = 1:size(n_sinais,1)
    % Numero e extens�o do arquivo
    arqnum =sprintf('%3d.mat',n_sinais(k0));
    load (arqnum); 
    N   = size(ecgs,1);     %Numero de amostras
    Ntypes=numel(type);     %Numero de tipos de anota��o

    %Plot 2D version of signal and labels
    figure; 
    maxEcg1 =max(ecgs(:,1))*ones(Ntypes,1);
    subplot(2,1,1); plot(ts(1:N),ecgs(1:N,1));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 1']); hold on;
                    plot(ts(ann(ann<N)+1),ecgs(ann(ann<N)+1),'ro');
                    text(ts(ann(ann<N)+1),maxEcg1,type);                    
    subplot(2,1,2); plot(ts(1:N),ecgs(1:N,2));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 2' ]);
    drawnow;

    Ts = 1/Fs;

    %% Tratamento de dados 

    [row_qrs,col_qrs] = find(type~='+');    %Posi��o das anota��es encontradas desconsiderar "non-beat annotations"
    qrs_ann = ann(row_qrs); %Vetor com as posi��es de instante das anota��es
    qrs_ann_t = qrs_ann*Ts; % Vetor com a poosi��o temporal das anota��es de pico R
    qrs_type = type(row_qrs,1); %Vetor com as etiquetas das anota��es

    % Caso seja necess�rio um pr�-processamento no sinal 
    % Normaliza��o no caso
    ecg = ecgs(1:N,1)/max(abs(ecgs(1:N,1)));

    %% Implementa��o do filtro passa-baixas

    b_l = zeros(1,13);
    b_l(1) = 1; b_l(7) = -2; b_l(13) = 1;
    a_l = 32*[1 -2 1];
    ecg_l = filter(b_l,a_l,ecg);

    % Butterworth: passa alta, ordem 8, fc = 5Hz e fa = 360Hz
    % fc = 11; % [Hz]
    % 
    % [b_l,a_l] = butter(12,fc/(Fs/2),'low'); % Utiliza a frequ�ncia de corte normalizada
    % ecg_l = filter(b_h,a_h,ecg);

    %% Implementa��o do filtro passa-altas

    % b_h = zeros(1,33);
    % b_h(1) = -1; b_h(17) = 32; b_h(18)= -32 ; b_h(33) = 1;
    % a_h = 32*[1 -1];
    % % Resultado filtrado em um passa banda (band-pass)
    % ecg_bp = filter(b_h,a_h,ecg_l);

    % Butterworth: passa alta, ordem 8, fc = 5Hz e fa = 360Hz
    % Defini��o das frequ�ncias
    fc = 2.75; % [Hz] 3 para o sinal 119

    [b_h,a_h] = butter(2,fc/(Fs/2),'high'); % Utiliza a frequ�ncia de corte normalizada (a ordem afeta muito o resultado)
    ecg_bp = filter(b_h,a_h,ecg_l);
    % Quanto menor a ordem, menor o atrado e menor a oscila��o no sinal filtrado

    % Retirou as oscila��es de baixiss�ma frequ�ncia, aproximando a m�dia de um n�vel DC nulo.

    %% FFT Filtros, OBSERVA��O DA FREQ DE CORTE, 5-11 Hz

    % obs.: O filtro foi originalmente projetado em uma fs = 200 Hz (diferente da utilizada aqui, aten��o na analise da resp. freq.)
    % Para mais detalher, ver o artigo corrigido pelo Tompkins
    fs = 200;

    % PASSA-BAIXA
    % Busca pela frequ�ncia de corte
    [h,w_Norm_l] = freqz(b_l,a_l,2048);
    mag = 20*log10(abs(h));
    w_piNorm_h = w_Norm_l/pi; % normalizada por pi
    for k1 = 1:size(mag,1)
        if (mag(k1) <= -3)
            break
        end
    end
    % Determina��o da frequencia do LP
    cutoff_freq_Norm = w_Norm_l(k1)/pi;
    cutoff_freq_l = cutoff_freq_Norm*fs/2;

    % PASSA-ALTA
    % Busca pela frequ�ncia de corte
    [h,w_Norm_h] = freqz(b_h,a_h,2048);
    mag = 20*log10(abs(h));
    w_piNorm_h = w_Norm_h/pi; % normalizada por pi
    for k1 = 1:size(mag,1)
        if (mag(k1) >= -3)
            break
        end
    end
    % Determina��o da frequencia do HP
    cutoff_freq_Norm_h = w_Norm_h(k1)/pi;
    cutoff_freq_h = cutoff_freq_Norm_h*Fs/2;

    %% Deriva��o

    ecg_d = zeros(size(ecg_bp));

    for k1 = 1:size(ecg_bp,1)
        if k1 > 4
            ecg_d(k1) = (1/8)*(2*ecg_bp(k1) + ecg_bp(k1-1) - ecg_bp(k1-3) - 2*ecg_bp(k1-4));
        elseif k1 == 3
            ecg_d(k1) = (1/8)*(2*ecg_bp(k1) + ecg_bp(k1-1));
        elseif k1 == 2
            ecg_d(k1) = (1/8)*(2*ecg_bp(k1));
        elseif k1 == 1
            ecg_d(k1) = (1/8)*(2*ecg_bp(k1));
        end
    end

    %% Squaring

    ecg_square = ecg_d.^2;

    %% Integra��o
    % M�todo utilidado: M�dia m�vel com janela igual a 30 elementos

    int_N = 30;
    int_mask = (1/int_N)*ones(int_N,1);

    % Convolu��o com os valores centrais, sem os resultantes do "zero padding".
    ecg_int = conv(ecg_square,int_mask,'same');

    %% Algoritmo de Pan-Tompkins

    % Definic�o da matriz de eventos
    % tamanho exagerado, reber� o valor 1 na posi��o vetorial equivalente ao instante de tempo em que ocorrer um evento. 
    eventos = NaN(size(ecg_int));

    % Sinal pr�-processado
    ecg_pp = ecg_int;

    % Vetores de m�dia
    RR_AVERAGE1_ARRAY = zeros(8,1);
    RR_AVERAGE2_ARRAY = zeros(8,1);

    % Medias
    RR_AVERAGE1 = 0;
    RR_AVERAGE2 = 0;

    % Limites
    RR_LOW_LIMIT = 0;
    RR_HIGH_LIMIT = 0;

    % Indicadores temporais tempor�rios de eventos 
    RR_temp1 = 0;
    RR_temp2 = 0;
    % obs.: RR_ind1 e RR_ind2 devem ir se alternando para 

    % Tamano do intervalo RR [s]
    RR_intervalo = 0;

    % Inicializando o limiares, utilizando o primeiro 1 segundo de dados (treinamento do algoritmo)
    N_treino = Fs; % N�mero de amostras em um segundo (tempo normalmente utilizado para treinar o algoritmo)

    % Pico de sinal
    SPKI = max(ecg_pp(1:N_treino));

    % Pico de ru�do
    NPKI = 0;

    % Inicializando o contador de pontos acima do threshold em um evento integrado
    % Lembrando que o �nico ponto que nos interessa � o disparo, o instante inicial do evento com ecg_pp(k1) > THRESHOLD_I1
    % A integra��o cria um evento mais longo que o pico R (ver Rangayyan pg. 190)
    cont_above_thr = 0;

    % N�mero de eventos total
    evento_cont = 0;

    THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
    THRESHOLD_I2 = 0.5*THRESHOLD_I1;

    THRESHOLD_I1_fn = NPKI + 0.25*(SPKI - NPKI);
    THRESHOLD_I2_fn = 0.5*THRESHOLD_I1;

    % Detec��o de eventos sem detec��o de falsos negativos
    % Aplica��o no sinal do tempo de treinamento em diante
    for k1 = (N_treino+1):size(ecg_pp,1)  

        if (ecg_pp(k1) > THRESHOLD_I1) 
            SPKI = 0.125*ecg_pp(k1) + 0.875*SPKI; 

            % Marca��o apenas do disparo do pico R
            if(cont_above_thr < 1) 
                eventos(k1) = 1;
            end
            cont_above_thr = cont_above_thr + 1;

            if(eventos(k1) == 1 && isnan(eventos(k1-1))) % Contador zerado a cada novo evento
                % Incremento de ocorr�ncia de evento
                evento_cont = evento_cont + 1;

                RR_temp1 = k1; % associando o instante de tempo inferior

                % C�lculo do intervalo RR
                RR_intervalo = (RR_temp1-RR_temp2)*Ts; % em segundos j�

                if (evento_cont > 1) % Para que se tenha realmente uma diferen�a temporal entre picos R-R 
                    % Guardando o tamanho do intervalo no vetor antes de zerar (m�todo de array push-remove)
                    RR_AVERAGE1_ARRAY(1:end-1) = RR_AVERAGE1_ARRAY(2:end);
                    RR_AVERAGE1_ARRAY(end) = RR_intervalo; % Associando a diferen�a de tempo
                    % C�lculo da m�dia considerando apenas os elementos n�o nulos
                    RR_AVERAGE1 = sum(RR_AVERAGE1_ARRAY)/nnz(RR_AVERAGE1_ARRAY);% mean(RR_AVERAGE1_ARRAY); % C�lculo da m�dia 1

                    % No inicio todos s�o nulos, ent�o deve haver um push inicial para esse caso
                    if (all(~RR_AVERAGE2_ARRAY) || (RR_AVERAGE1 > RR_LOW_LIMIT && RR_AVERAGE1 < RR_HIGH_LIMIT))
                        RR_AVERAGE2_ARRAY(1:end-1) = RR_AVERAGE2_ARRAY(2:end);
                        RR_AVERAGE2_ARRAY(end) = RR_AVERAGE1; % RR_AVERAGE2_ARRAY recebe o mesmo intervalo RR calculado para RR_AVERAGE1
                        % C�lculo da m�dia considerando apenas os elementos n�o nulos
                        RR_AVERAGE2 = sum(RR_AVERAGE2_ARRAY)/nnz(RR_AVERAGE2_ARRAY); % C�lculo da m�dia 2
                        % Recalculando os limites 
                        RR_LOW_LIMIT = 0.92*RR_AVERAGE2;
                        RR_HIGH_LIMIT = 1.16*RR_AVERAGE2;
                    end
                end           
            end       
        else  
            NPKI = 0.125*ecg_pp(k1) + 0.875*NPKI;
            % Significa que acabou a zona acima do thr gerada pela integra��o, logo deve-se zerar o contador de anota��es acima do thr
            cont_above_thr = 0; 
        end

        % Contador de anota��es acima do thr geradas pela integra��o
        cont_above_thr_sb = 0;
        % Searchback
        % Checando se o intervalo RR calculado na itera��o atual ultrapassa o limite RR_MISSED_LIMIT = 166%RR_AVERAGE2
        if (RR_intervalo > 1.66*RR_AVERAGE2 && evento_cont > 1) 
            % A segunda parte da condicional indica que agora sim trata-se de uma diferen�a RR
            fprintf('searchback \n');      
            for k2 = RR_temp2:k1
                if(ecg_pp(k2) > THRESHOLD_I2 && cont_above_thr_sb < 1) 
                    % A segunda condicional garante que s� o disparo do pico R vai ser anotado
                    % Instante recebe dete��o de evento baseado em THRESHOLD_I2, mas nada � atualizado pois esse � um outlier
                    eventos(k2) = 1;
                    cont_above_thr_sb = cont_above_thr_sb + 1;
                end
            end
            cont_above_thr_sb = 0;
        end

        % Atualiza��o do tempor�rio 2 (atrasada par autiliza��o no if anterior) 
        RR_temp2 = RR_temp1;

        % Atualizar os thresholds
        THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
        THRESHOLD_I2 = 0.5*THRESHOLD_I1;

    end

    %% Dete��o de falsos eventos 
    % Ser�o comparados os eventos com uma janela de toler�ncia igual a 50ms
    % ref.: The Accuracy on the Common Pan-Tompkins Based QRS Detection Methods Through Low-Quality Electrocardiogram Database
    % Chengyu Liu et. al.
    % obs.: Cabe lembrar que a compara��o deve descartar o primeiro segundo da s�ria, por conta do treino do algoritmo.
    % obs2.: O Butterworth implementado com a fun��o butter insere um atraso maior logo, usar 100ms de cada lado da janela
    tolerancia = 100e-3; % [s]

    % Cria��o de um vetor com a posi��o temporal dos eventos (baseado no vetor eventos)
    eventos_tpos  = zeros(evento_cont,1);
    % Inicializando o segundo contador para atribui��o do instante
    k2 = 1;
    for k1 = 1:size(eventos,1)
        if ~isnan(eventos(k1)) % Toda vez que o elemento n�o for NaN
            eventos_tpos(k2) = k1*Ts;
            k2 = k2 + 1;
        end
    end

    % Cria��o do vetor de eventos sem o primeiro segundo. (em segundos, baseado nas anota��es)
    qrs_comparativo = qrs_ann_t(qrs_ann_t > 1);

    % Fazer a janela com uma matriz de duas colunas e n�mero de eventos linhas (contendo os extremos de toler�ncia)
    janela_tol = zeros(size(qrs_comparativo,1),2);
    for k1 = 1:size(qrs_comparativo,1)
        janela_tol(k1,1) = qrs_comparativo(k1) - tolerancia; % Extremo inferior de  toler�ncia para o evento 
        janela_tol(k1,2) = qrs_comparativo(k1) + tolerancia; % Extremo superior de  toler�ncia para o evento 
    end

    % Vetor para determinar os falsos positivos 
    % NaN: n�o h� intervalo associado ao evento
    % 1: h� intervalo associado ao evento
    fp_eventos = NaN(size(eventos_tpos,1),1);

    % Vetor para determinar os falsos negativos 
    % NaN: n�o h� evento associado ao intervalo
    % 1: h� evento associado ao intervalo
    fn_eventos = NaN(size(janela_tol,1),1);

    % Falsos positivos (checar se cada evento est� abrigado em um intervalo)
    for k1 = 1:size(eventos_tpos,1)
        for k2 = 1:size(janela_tol,1) 
            if (eventos_tpos(k1) > janela_tol(k2,1) && eventos_tpos(k1) < janela_tol(k2,2))
                % Associar evento ao seu intervalo
                fp_eventos(k1) = 1;
            end
        end
    end
    % C�lculo da quantidade de falsos positivos 
    FP = size(fp_eventos,1) - sum(fp_eventos,'omitnan'); % False positive
    % C�lculo da quantidade de verdadeiros positivos
    TP = size(fp_eventos,1) - FP; % True positive

    % Falsos negativos (checar se cada intervalo abriga um evento)
    for k1 = 1:size(janela_tol,1) 
        for k2 = 1:size(eventos_tpos,1)
            if (eventos_tpos(k2) > janela_tol(k1,1) && eventos_tpos(k2) < janela_tol(k1,2))
                % Associar evento ao seu intervalo
                fn_eventos(k1) = 1;
            end
        end
    end
    % C�lculo da quantidade de falsos negativos 
    FN = size(fn_eventos,1) - sum(fn_eventos,'omitnan'); % False negative
    % C�lculo da quantidade de verdadeiros negativos 
    TN = size(fn_eventos,1) - FN; % True negative
    %% C�lculo das m�tricas
    % Colocar as m�tricas em uma table e associar os valores calculados inicialmente a vetores, uma vez pra cada sinal

    %Sinal(k0) = n_sinais(k0);
    Media(k0) = mean(ecgs(1:N,1));
    DP(k0) = std(ecgs(1:N,1));
    Media_norm(k0) = mean(ecg);
    DP_norm(k0) = std(ecgs(1:N,1));
    Nv(k0) = size(qrs_comparativo,1);
    Nd(k0) = evento_cont;
    FP_T(k0) = FP;
    p_FP(k0) = (FP/Nv(k0))*100;
    FN_T(k0) = FN;
    p_FN(k0) = (FN/Nv(k0))*100;

    % Relacionados ao IRR
    % Cria��o de um vetor com os intervalor RR
    IRR = diff(qrs_comparativo);
    Media_IRR(k0) = mean(IRR);
    DP_IRR(k0) = std(IRR);
    
end

%% Montagem da tabela
% Associa��o dos vetores as posi��es da tabela 
T_metricas = table(Sinal, Media, DP, Media_norm, DP_norm, Nv, Nd, FP_T, p_FP, FN_T, p_FN, Media_IRR, DP_IRR);

% Salvar em xlsx para facilitar a convers�o para LaTex table
%writetable(T_metricas,'tabela_metricas.xlsx');




