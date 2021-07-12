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

%%
for n =111:111
    arqnum =sprintf('%3d.mat',n);
    load (arqnum);      
       % save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as vari�veis
       % ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna cont�m sinal com Nt amostras j� calibradas em mV
       % ts(Nt,1): [double] vetor com Nt instantes de tempo em s
       % Fs: frequencia de amostragem em Hz
       % ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
       %       anotacao ocorreu na amostra 270, entao p=ann(3) retornar� p=270.
       % type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.      
    N   =size(ecgs,1);
    fprintf('\n %s com %d amostras',arqnum,N);
    Ntypes=numel(type);
    %Plot 2D version of signal and labels
    figure; 
    maxEcg1 =max(ecgs(:,1))*ones(Ntypes,1);
    subplot(2,1,1); plot(ts(1:N),ecgs(1:N,1));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 1']); hold on;
                    plot(ts(ann(ann<N)+1),ecgs(ann(ann<N)+1),'ro');
                    text(ts(ann(ann<N)+1),maxEcg1,type);                    
    subplot(2,1,2); plot(ts(1:N),ecgs(1:N,2));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 2' ]);
    drawnow;
end
fprintf('\n Fim \n');
Ts = 1/Fs;
%% Load do sinal no espa�o de trabalho

n =118;
% Numero e extens�o do arquivo
arqnum =sprintf('%3d.mat',n);
load (arqnum); 
N   = size(ecgs,1);     %Numero de amostras
Ntypes=numel(type);     %Numero de tipos de anota��o

%% Plot do sinal

figure; 
plot(ts(1:N),ecgs(1:N,1));
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1');    
drawnow;

%% Tratamento de dados 

%type = ['A';'B';'V';'N';'A';'B'] % TESTE
%ann = [1;2;3;4;5;6]% TESTE
[row_qrs,col_qrs] = find(type~='N' & type~='V');    %Posi��o das anota��es encontradas
qrs_ann = ann(row_qrs);     %Vetor com as posi��es das anota��es
qrs_type = type(row_qrs,1); %Vetor com as etiquetas das anota��es

ecg = ecgs(1:N,1) - mean(ecgs(1:N,1)); % N�O ENTENDI (ROBERTO)

%% Implementa��o do filtro passa-baixas

b_l = zeros(1,13);
b_l(1) = 1; b_l(7) = -2; b_l(13) = 1;
a_l = 32*[1 -2 1];
ecg_l = filter(b_l,a_l,ecg);

%% Implementa��o do filtro passa-altas

b_h = zeros(1,33);
b_h(1) = -1; b_h(17) = 32; b_h(18)= -32 ; b_h(33) = 1;
a_h = 32*[1 -1];
% Resultado filtrado em um passa banda (band-pass)
ecg_bp = filter(b_h,a_h,ecg_l);

% Retirou as oscila��es de baixiss�ma frequ�ncia, aproximando a m�dia de um n�vel DC nulo.

%% FFT Filtros, OBSERVA��O DA FREQ DE CORTE, 5-11 Hz
figure;
freqz(b_l,a_l)
title('Resposta em frequ�ncia do filtro passa-baixas')

figure;
freqz(b_h,a_h)
title('Resposta em frequ�ncia do filtro passa-altas')

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
eventos_fn = NaN(size(ecg_int));

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
        eventos(k1) = 1;
        cont_above_thr = cont_above_thr + 1;
        
        if(eventos(k1) == 1 && isnan(eventos(k1-1))) % Contador zerado a cada novo evento
            
            % Incremento de ocorr�ncia de evento
            evento_cont = evento_cont + 1;
            
            RR_temp1 = k1*Ts; % associando o instante de tempo inferior

            % C�lculo do intervalo RR
            RR_intervalo = RR_temp1-RR_temp2;
            
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
        eventos(k1-(cont_above_thr-1):k1) = NaN; % Atribuindo NaN a todos os instantes do evento integrado, menos o primeiro (disparo)
        cont_above_thr = 0;
        
        % Atualiza��o do indicador tempor�rio
        RR_temp2 = RR_temp1;
    end
          
    % Checando se o intervalo RR calculado na itera��o atual ultrapassa o limite RR_MISSED_LIMIT = 166%RR_AVERAGE2
    if (RR_intervalo > 1.66*RR_AVERAGE2 && evento_cont > 1) 
        % A segunda parte da condicional indica que agora sim trata-se de uma diferen�a RR
        % IMPLEMENTAR o searchback
        fprintf('searchback \n');
    end
        
        
    % Atualizar os thresholds
    THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
    THRESHOLD_I2 = 0.5*THRESHOLD_I1;
    
end


%% Plot comparativo com sele��o do per�odo QRS

ann_number = 100;
t1 = ann(ann_number,1);
t2 = ann(ann_number+3,1);

figure; 
subplot(5,1,1); 
plot(ts(t1:t2),ecg(t1:t2));
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal ECG original'); 
hold on;
subplot(5,1,2); 
plot(ts(t1:t2),ecg_bp(t1:t2));
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal ECG filtrado');
subplot(5,1,3); 
plot(ts(t1:t2),ecg_d(t1:t2));
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal ECG derivado'); 
subplot(5,1,4); 
plot(ts(t1:t2),ecg_square(t1:t2));
ylabel('mV^2'); 
xlabel('s'); 
grid on;
title('Sinal ECG squared'); 
subplot(5,1,5); 
plot(ts(t1:t2),ecg_int(t1:t2));
ylabel('mV^2'); 
xlabel('s'); 
grid on;
title('Sinal ECG integrado'); 
drawnow;


%% Plots 

% Plot filtragem
figure; 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_bp);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado');    
drawnow;

% Deriva��o
figure() 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_d);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado e derivado');    
drawnow;

% Squared
figure() 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_square);
ylabel('mV^2'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado, derivado e squared');    
drawnow;

% Integrado
figure() 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_int);
ylabel('mV^2'); 
xlabel('s'); 
grid on;
title('Sinal 1 filtrado, derivado, squared e integrado');    
drawnow;

%%
% Detec��o
figure; 
plot(ts(1:N),ecg,ts(1:N),eventos,'r*');
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 e a dete��o de eventos');    
drawnow;
% obs.: deve-se pegar apenas o primeiro 1 para cada evento, pois os outros s�o parte do pronlogamento do sinal dada a integra��o

%%
ann_number =40;
t1 = ann(ann_number,1);
t2 = ann(ann_number+20,1);
qrs = ecgs(t1:t2,1);
qrs_t = ts(t1:t2);

figure; 
subplot(2,1,1)
plot(ts(t1:t2),ecg(t1:t2),ts(t1:t2),eventos(t1:t2),'r*');
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 e a dete��o de eventos');  
subplot(2,1,2)
plot(ts(t1:t2),ecg(t1:t2),ts(t1:t2),eventos_fn(t1:t2),'r*');
ylabel('mV'); 
xlabel('s'); 
grid on;
title('Sinal 1 e a dete��o de eventos com detec��o de falsos negativos'); 
drawnow;

figure; 
plot(eventos==eventos_fn)
title('Falsos negativos'); 
