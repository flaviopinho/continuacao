%inicializacao.m
function[u0, t0, estbilidade0, h0, w0, modelo, parametros, id_saida]=inicializacao()

clc;
close('all');
fclose('all');

% Numero de variaveis
modelo.n=2;
% Numero de parametros
modelo.p=3;
% Funcao residuo
modelo.residuo=@residuo_predador_presa;
% Funcao jacobiana
modelo.jacobiana=@jacobiana_predador_presa;
% Funcao de determinacao da estabilidade
modelo.estabilidade=@estabilidade_ponto_fixo;
% Limites de interesse das variaveis
modelo.limites=[ -1,   2;
                 -1,   3;
                2.9, 3.1;
                4.9, 5.1;
                  0,   1];

% Determinacao de um ponto regular inicial
u0=[1/3; 2; 3; 5; 0];
H0=residuo_predador_presa(u0);
J0=jacobiana_predador_presa(u0);
t0=vetor_tangente(J0);
estbilidade0=estabilidade_ponto_fixo(u0, J0, modelo);

% Passo e sentido iniciais
h0=0.0001;
w0=1;
       
% Numero de pontos maximo determinados
parametros.cont_max=1000;
% Numero maximo de iteracoes do corretor
parametros.i_max=6;
% Tamanho de passo minimo e maximo
parametros.h_min=1e-6;
parametros.h_max=0.05;
% Parametro de ajuste do passo
parametros.i_controle=3;
% Tolerancias
parametros.tol1=1e-6;
parametros.tol2=1e-6;

% Identificador do arquivo de saida de dados
id_saida=fopen('resultados_predador_presa1.txt', 'w');

end