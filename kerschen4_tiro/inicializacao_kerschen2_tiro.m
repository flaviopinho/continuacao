function [u0, t0, estbilidade0, h0, w0, modelo, parametros, id_saida]=inicializacao_kerschen2_tiro()

%% Inicialização do balanço harmônico

addpath('Tiro múltiplo');
tiro=tiro_multiplo;
tiro.din_func=@(x, pars) kerschen4(1, x, pars);
tiro.din_jacobiana_x=@(x, pars) jacobiana_kerschen4_x(1, x, pars);
tiro.din_jacobiana_pars=@(x, pars) jacobiana_kerschen4_pars(1, x, pars);
        
tiro.autonomo=true;
tiro.complexo=false;
        
tiro.n_din=4;
tiro.p=3;

tiro.Nt=10;

tiro=tiro.montar_modelo();

%% 
% Numero de variaveis
modelo.n=tiro.n;
% Numero de parametros
modelo.p=tiro.p;
% Funcao residuo
modelo.residuo=@(u) residuo_kerschen4_tiro_autonomo(tiro, u);
% Funcao jacobiana
modelo.jacobiana=@(u) jacobiana_kerschen4_tiro_autonomo(tiro, u);
% Funcao de determinacao da estabilidade
modelo.estabilidade=@estabilidade_tiro_Floquet;
% Limites de interesse das variaveis             
limites_Q=repmat([-Inf, Inf], tiro.n, 1);
        
limites_pars=[ -10,   10;
               -10,   10;
                 0,  0.45*2*pi;];
modelo.limites=[limites_Q; limites_pars];
modelo.outdof1=modelo.n+modelo.p;
modelo.outdof2=5;
% pasta do modelo
modelo.folder='kerschen4_tiro';
modelo.data=tiro;

%% Parametros iniciais
x0=ones(tiro.n_din, 1)*0.001;
pars0=[0; 0.00001; 1.001];
X0=ponto_inicial_tiro(tiro, x0, pars0);

u0=[X0; pars0];

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
id_saida=fopen(fullfile(modelo.folder, 'resultados_kerschen2_tiro_1.txt'), 'w');

save(fullfile(modelo.folder, 'kerschen2_tiro.mat'), 'modelo', 'parametros');

end
