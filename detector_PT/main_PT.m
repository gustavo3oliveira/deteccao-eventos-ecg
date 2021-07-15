% Author: Gustavo Henrique de Oliveira, Maria Antònia Marzo Barceló, Roberto Costa Ceccato
% Date: 27/06/2021
% Description: Script principal para o projeto da Parte B (Prof. Sérgio)

% save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as variáveis
% ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna contém sinal com Nt amostras já calibradas em mV
% ts(Nt,1): [double] vetor com Nt instantes de tempo em s
% Fs: frequencia de amostragem em Hz
% ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
% anotacao ocorreu na amostra 270, entao p=ann(3) retornará p=270.
% type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.      

close all 

%% Load do sinal no espaço de trabalho
% save(arq_mat,'ecgs','ts','Fs','ann','type');  Ver documentacao abaixo sobre as variáveis
% ecgs(Nt,2):[double] matriz com Nt linhas por 2 colunas. Cada coluna contém sinal com Nt amostras já calibradas em mV
% ts(Nt,1): [double] vetor com Nt instantes de tempo em s
% Fs: frequencia de amostragem em Hz
% ann(Na,1):[double] vetor contendo a posicao (em amostra) da anotacao. Tem 'Na' anotacoes. Ex.: suponha que a terceira 
%       anotacao ocorreu na amostra 270, entao p=ann(3) retornará p=270.
% type(Na,1):[char] vetor com a anotacao para cada uma das 'Na' anotacoes.  

n =119;
% Numero e extensão do arquivo
arqnum =sprintf('%3d.mat',n);
load (arqnum); 
N   = size(ecgs,1);     %Numero de amostras
Ntypes=numel(type);     %Numero de tipos de anotação

%Plot 2D version of signal and labels
figure; 
maxEcg1 =max(ecgs(:,2))*ones(Ntypes,1);
subplot(2,1,1); plot(ts(1:N),ecgs(1:N,2));

                ylabel('mV'); xlabel('s'); 
                grid on;
                title([arqnum ' sinal 1']); 
                hold on;
                
                plot(ts(ann(ann<N)+1),ecgs(ann(ann<N)+1),'ro');
                text(ts(ann(ann<N)+1),maxEcg1,type);                    
subplot(2,1,2); plot(ts(1:N),ecgs(1:N,2));ylabel('mV'); xlabel('s'); grid on;title([arqnum ' sinal 2' ]);
drawnow;

Ts = 1/Fs;




%% Tratamento de dados 

[row_qrs,col_qrs] = find(type~='+');    %Posição das anotações encontradas desconsiderar "non-beat annotations"
qrs_ann = ann(row_qrs); %Vetor com as posições de instante das anotações
qrs_ann_t = qrs_ann*Ts; % Vetor com a poosição temporal das anotações de pico R
qrs_type = type(row_qrs,1); %Vetor com as etiquetas das anotações

% Caso seja necessário um pré-processamento no sinal 
% Normalização no caso
ecg = ecgs(1:N,2)/max(abs(ecgs(1:N,1)));

%% Plot do sinal

f = figure;
f.Position = [100 100 900 400]; 
subplot(2,1,1); 
plot(ts(1:N),ecgs(1:N,2));
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)');
title('Sinal de electrocardiograma original');    
grid on;
subplot(2,1,2); 
drawnow;
plot(ts(1:N),ecg(1:N));
ylabel('Voltagem (mV)');
ylim([-3 3])
xlabel('Tempo (s)');
title('Sinal de electrocardiograma normalizado'); 
grid on;
drawnow;
saveas(f,'Sinal 119 original e normalizado.png')


%% Implementação do filtro passa-baixas

b_l = zeros(1,13);
b_l(1) = 1; b_l(7) = -2; b_l(13) = 1;
a_l = 32*[1 -2 1];
ecg_l = filter(b_l,a_l,ecg);

%% Implementação do filtro passa-altas

b_h = zeros(1,33);
b_h(1) = -1; b_h(17) = 32; b_h(18)= -32 ; b_h(33) = 1;
a_h = 32*[1 -1];
% Resultado filtrado em um passa banda (band-pass)
ecg_bp = filter(b_h,a_h,ecg_l);

% Retirou as oscilações de baixissíma frequência, aproximando a média de um nível DC nulo.

%% FFT Filtros, OBSERVAÇÃO DA FREQ DE CORTE, 5-11 Hz
f = figure;
f.Position = [100 100 800 300]; 
freqz(b_l,a_l)
title('Resposta em frequência do filtro passa-baixas')
saveas(f,'Filtro passa baixas original.png')

f = figure;
f.Position = [100 100 800 300]; 
freqz(b_h,a_h)
title('Resposta em frequência do filtro passa-altas')
saveas(f,'Filtro passa altas original.png')


%% Derivação

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

%% Integração
% Método utilidado: Média móvel com janela igual a 30 elementos

int_N = 30;
int_mask = (1/int_N)*ones(int_N,1);

% Convolução com os valores centrais, sem os resultantes do "zero padding".
ecg_int = conv(ecg_square,int_mask,'same');

%% Algoritmo de Pan-Tompkins

% Definicão da matriz de eventos
% tamanho exagerado, reberá o valor 1 na posição vetorial equivalente ao instante de tempo em que ocorrer um evento. 
eventos = NaN(size(ecg_int));

% Sinal pré-processado
ecg_pp = ecg_int;

% Vetores de média
RR_AVERAGE1_ARRAY = zeros(8,1);
RR_AVERAGE2_ARRAY = zeros(8,1);

% Medias
RR_AVERAGE1 = 0;
RR_AVERAGE2 = 0;

% Limites
RR_LOW_LIMIT = 0;
RR_HIGH_LIMIT = 0;

% Indicadores temporais temporários de eventos 
RR_temp1 = 0;
RR_temp2 = 0;
% obs.: RR_ind1 e RR_ind2 devem ir se alternando para 

% Tamano do intervalo RR [s]
RR_intervalo = 0;

% Inicializando o limiares, utilizando o primeiro 1 segundo de dados (treinamento do algoritmo)
N_treino = Fs; % Número de amostras em um segundo (tempo normalmente utilizado para treinar o algoritmo)

% Pico de sinal
SPKI = max(ecg_pp(1:N_treino));

% Pico de ruído
NPKI = 0;

% Inicializando o contador de pontos acima do threshold em um evento integrado
% Lembrando que o único ponto que nos interessa é o disparo, o instante inicial do evento com ecg_pp(k1) > THRESHOLD_I1
% A integração cria um evento mais longo que o pico R (ver Rangayyan pg. 190)
cont_above_thr = 0;

% Número de eventos total
evento_cont = 0;

THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
THRESHOLD_I2 = 0.5*THRESHOLD_I1;

THRESHOLD_I1_fn = NPKI + 0.25*(SPKI - NPKI);
THRESHOLD_I2_fn = 0.5*THRESHOLD_I1;

% Detecção de eventos sem detecção de falsos negativos
% Aplicação no sinal do tempo de treinamento em diante
for k1 = (N_treino+1):size(ecg_pp,1)  
  
    if (ecg_pp(k1) > THRESHOLD_I1) 
        SPKI = 0.125*ecg_pp(k1) + 0.875*SPKI; 
        eventos(k1) = 1;
        cont_above_thr = cont_above_thr + 1;
        
        if(eventos(k1) == 1 && isnan(eventos(k1-1))) % Contador zerado a cada novo evento
            
            % Incremento de ocorrência de evento
            evento_cont = evento_cont + 1;
            
            RR_temp1 = k1; % associando o instante de tempo inferior

            % Cálculo do intervalo RR
            RR_intervalo = (RR_temp1-RR_temp2)*Ts; % em segundos já
            
           
            if (evento_cont > 1) % Para que se tenha realmente uma diferença temporal entre picos R-R 
                % Guardando o tamanho do intervalo no vetor antes de zerar (método de array push-remove)
                RR_AVERAGE1_ARRAY(1:end-1) = RR_AVERAGE1_ARRAY(2:end);
                RR_AVERAGE1_ARRAY(end) = RR_intervalo; % Associando a diferença de tempo
                % Cálculo da média considerando apenas os elementos não nulos
                RR_AVERAGE1 = sum(RR_AVERAGE1_ARRAY)/nnz(RR_AVERAGE1_ARRAY);% mean(RR_AVERAGE1_ARRAY); % Cálculo da média 1


                % No inicio todos são nulos, então deve haver um push inicial para esse caso
                if (all(~RR_AVERAGE2_ARRAY) || (RR_AVERAGE1 > RR_LOW_LIMIT && RR_AVERAGE1 < RR_HIGH_LIMIT))
                    RR_AVERAGE2_ARRAY(1:end-1) = RR_AVERAGE2_ARRAY(2:end);
                    RR_AVERAGE2_ARRAY(end) = RR_AVERAGE1; % RR_AVERAGE2_ARRAY recebe o mesmo intervalo RR calculado para RR_AVERAGE1
                    % Cálculo da média considerando apenas os elementos não nulos
                    RR_AVERAGE2 = sum(RR_AVERAGE2_ARRAY)/nnz(RR_AVERAGE2_ARRAY); % Cálculo da média 2
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
        
    end
    
    cont_above_thr_sb = 0;
    % Searchback
    % Checando se o intervalo RR calculado na iteração atual ultrapassa o limite RR_MISSED_LIMIT = 166%RR_AVERAGE2
    if (RR_intervalo > 1.66*RR_AVERAGE2 && evento_cont > 1) 
        % A segunda parte da condicional indica que agora sim trata-se de uma diferença RR
        fprintf('searchback \n');      
        for k2 = RR_temp2:k1
            if(ecg_pp(k2) > THRESHOLD_I2 && cont_above_thr_sb < 1) 
                % A segunda condicional garante que só o disparo do pico R vai ser anotado
                % Instante recebe deteção de evento baseado em THRESHOLD_I2, mas nada é atualizado pois esse é um outlier
                eventos(k2) = 1;
                cont_above_thr_sb = cont_above_thr_sb + 1;
            end
        end
        cont_above_thr_sb = 0;
    end

    % Atualização do temporário 2 (atrasada par autilização no if anterior) 
    RR_temp2 = RR_temp1;
   
    % Atualizar os thresholds
    THRESHOLD_I1 = NPKI + 0.25*(SPKI - NPKI);
    THRESHOLD_I2 = 0.5*THRESHOLD_I1;
    
end

%% Deteção de falsos eventos 
% Serão comparados os eventos com uma janela de tolerância igual a 50ms
% ref.: The Accuracy on the Common Pan-Tompkins Based QRS Detection Methods Through Low-Quality Electrocardiogram Database
% Chengyu Liu et. al.
% obs.: Cabe lembrar que a comparação deve descartar o primeiro segundo da séria, por conta do treino do algoritmo.
tolerancia = 50e-3; % [s]

% Criação de um vetor com a posição temporal dos eventos (baseado no vetor eventos)
eventos_tpos  = zeros(evento_cont,1);
% Inicializando o segundo contador para atribuição do instante
k2 = 1;
for k1 = 1:size(eventos,1)
    if ~isnan(eventos(k1)) % Toda vez que o elemento não for NaN
        eventos_tpos(k2) = k1*Ts;
        k2 = k2 + 1;
    end
end

% Criação do vetor de eventos sem o primeiro segundo. (em segundos, baseado nas anotações)
qrs_comparativo = qrs_ann_t(qrs_ann_t > 1);

% Fazer a janela com uma matriz de duas colunas e número de eventos linhas (contendo os extremos de tolerância)
janela_tol = zeros(size(qrs_comparativo,1),2);
for k1 = 1:size(qrs_comparativo,1)
    janela_tol(k1,1) = qrs_comparativo(k1) - tolerancia; % Extremo inferior de  tolerância para o evento 
    janela_tol(k1,2) = qrs_comparativo(k1) + tolerancia; % Extremo superior de  tolerância para o evento 
end

% Vetor para determinar os falsos positivos 
% NaN: não há intervalo associado ao evento
% 1: há intervalo associado ao evento
fp_eventos = NaN(size(eventos_tpos,1),1);

% Vetor para determinar os falsos negativos 
% NaN: não há evento associado ao intervalo
% 1: há evento associado ao intervalo
fn_eventos = NaN(size(janela_tol,1),1);

% Falsos positivos (checar se cada evento está abrigado em um intervalo)
for k1 = 1:size(eventos_tpos,1)
    for k2 = 1:size(janela_tol,1) 
        if (eventos_tpos(k1) > janela_tol(k2,1) && eventos_tpos(k1) < janela_tol(k2,2))
            % Associar evento ao seu intervalo
            fp_eventos(k1) = 1;
        end
    end
end
% Cálculo da quantidade de falsos positivos 
FP = size(fp_eventos,1) - sum(fp_eventos,'omitnan'); % False positive
% Cálculo da quantidade de verdadeiros positivos
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
% Cálculo da quantidade de falsos negativos 
FN = size(fn_eventos,1) - sum(fn_eventos,'omitnan'); % False negative
% Cálculo da quantidade de verdadeiros negativos 
TN = size(fn_eventos,1) - FN; % True negative

%% Cálculo das métricas
% Colocar as métricas em uma table e associar os valores calculados inicialmente a vetores, uma vez pra cada sinal

Sinal = n;
Media = mean(ecgs(1:N,1));
DP = std(ecgs(1:N,1));
Media_norm = mean(ecg);
DP_norm = std(ecgs(1:N,1));
Nv = size(qrs_comparativo,1);
Nd = evento_cont;
FP_T = FP;
p_FP = (FP/Nv)*100;
FN_T = FN;
p_FN = (FN/Nv)*100;

% Relacionados ao IRR
% Criação de um vetor com os intervalor RR
IRR = diff(qrs_comparativo);
Media_IRR = mean(IRR);
DP_IRR = std(IRR);

T_metricas = table(Sinal, Media, DP, Media_norm, DP_norm, Nv, Nd, FP_T, p_FP, FN_T, p_FN, Media_IRR, DP_IRR);

%% Plot comparativo com seleção do período QRS

ann_number = 40;
t1 = ann(ann_number,1);
t2 = ann(ann_number+3,1);

f = figure;
f.Position = [100 100 900 600]; 
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
saveas(f,'Sinal 119 Diagrama comparativo.png')



%%
% Detecção
f = figure;
f.Position = [100 100 900 300]; 
plot(ts(1:N),ecg,ts(1:N),eventos,'r*');
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)');  
grid on;
title('Sinal 1 e a deteção de eventos');    
drawnow;
saveas(f,'Sinal 119 detecao eventos.png')

% obs.: deve-se pegar apenas o primeiro 1 para cada evento, pois os outros são parte do pronlogamento do sinal dada a integração

%%
ann_number =30;
t1 = ann(ann_number,1);
t2 = ann(ann_number+20,1);
qrs = ecgs(t1:t2,1);
qrs_t = ts(t1:t2);

f = figure;
f.Position = [100 100 800 350]; 
plot(ts(t1:t2),ecg(t1:t2),ts(t1:t2),eventos(t1:t2),'r*');
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Deteção de anotações para 20 eventos');  
saveas(f,'Sinal 119 detecao 20 eventos.png')


%% Plots construção

% Plot filtragem
f = figure;
f.Position = [100 100 800 350]; 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG normalizado sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_bp);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG normalizado e filtrado');    
drawnow;
saveas(f,'Sinal 119 fitrado.png')


% Derivação
f = figure;
f.Position = [100 100 800 350]; subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)');  
grid on;
title('Sinal de ECG normalizado sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_d);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG normalizado, filtrado e derivado');    
drawnow;
saveas(f,'Sinal 119 derivado.png')


% Squared
f = figure;
f.Position = [100 100 800 350]; 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG normalizado sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_square);
ylabel('mV^2'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG filtrado, derivado e squared');    
drawnow;
saveas(f,'Sinal 119 squared.png')

% Integrado
f = figure;
f.Position = [100 100 800 350]; 
subplot(2,1,1); 
plot(ts(1:N),ecg);
ylabel('Voltagem (mV)'); 
xlabel('Tempo (s)'); 
grid on;
title('Sinal de ECG normalizado sem filtrar'); 
hold on;
subplot(2,1,2); 
plot(ts(1:N),ecg_int);
ylabel('mV^2'); 
xlabel('Tempo (s)');  
grid on;
title('Sinal de ECG filtrado, derivado, squared e integrado');    
drawnow;
saveas(f,'Sinal 119 integrado.png')

