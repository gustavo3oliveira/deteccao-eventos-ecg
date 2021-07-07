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
Fa          =360;  %Hz
for n =118:118
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
%% Load do sinal no espa�o de trabalho

n =111;
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

ecg = ecgs(1:N,1) - mean(ecgs(1:N,1));

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
int_mask = (1/30)*ones(int_N,1);

% Convolu��o com os valores centrais, sem os resultantes do "zero padding".
ecg_int = conv(ecg_square,int_mask,'same');

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

