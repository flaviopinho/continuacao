function [u0, t0, estbilidade0, h0, w0, modelo, parametros, id_saida]=inicializacao_kerschen2()

%% Inicialização do balanço harmônico

addpath('Balanço Harmônico');
bh=balanco_harmonico_AFT;
bh.din_func=@(x, pars) kerschen4(1, x, pars);
bh.din_jacobiana_x=@(x, pars) jacobiana_kerschen4_x(1, x, pars);
bh.din_jacobiana_pars=@(x, pars) jacobiana_kerschen4_pars(1, x, pars);
        
bh.autonomo=true;
bh.complexo=false;
        
bh.n_din=4;
bh.p=3;

bh.H=15;
bh.Nt=2^10;

bh=bh.montar_modelo();

%% 
% Numero de variaveis
modelo.n=bh.n;
% Numero de parametros
modelo.p=bh.p;
% Funcao residuo
modelo.residuo=@(u) residuo_kerschen4_bh_autonomo(bh, u);
% Funcao jacobiana
modelo.jacobiana=@(u) jacobiana_kerschen4_bh_autonomo(bh, u);
% Funcao de determinacao da estabilidade
modelo.estabilidade=@estabilidade_bh_Hill;
% Limites de interesse das variaveis             
limites_Q=repmat([-Inf, Inf], bh.n, 1);
        
limites_pars=[ -10,   10;
               -10,   10;
                 0,  0.45*2*pi;];
modelo.limites=[limites_Q; limites_pars];
modelo.outdof1=modelo.n+modelo.p;
modelo.outdof2=5;
% pasta do modelo
modelo.folder='kerschen4';
modelo.data=bh;

%% Parametros iniciais
x0=ones(bh.n_din, 1)*0.001;
pars0=[0; 0.00001; 1.001];
Q0=ponto_inicial_bh(bh, x0, pars0);

u0=[Q0; pars0];

H0=modelo.residuo(u0);
J0=modelo.jacobiana(u0);
t0=vetor_tangente(J0);
estbilidade0=modelo.estabilidade(u0, J0, modelo);

%% parametros

% Passo e sentido iniciais
h0=0.01;
w0=1;
       
% Numero de pontos maximo determinados
parametros.cont_max=10000;
% Numero maximo de iteracoes do corretor
parametros.i_max=5;
% Tamanho de passo minimo e maximo
parametros.h_min=1e-6;
parametros.h_max=20;
% Parametro de ajuste do passo
parametros.i_controle=3;
% Tolerancias
parametros.tol1=1e-10;
parametros.tol2=1e-10;

% Identificador do arquivo de saida de dados
id_saida=fopen(fullfile(modelo.folder, 'resultados_kerschen2_1.txt'), 'w');

save(fullfile(modelo.folder, 'kerschen2.mat'), 'modelo', 'parametros');

end