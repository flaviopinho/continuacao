function [u0, t0, estbilidade0, h0, w0, modelo, parametros, id_saida]=inicializacao_duffing_tiro()

%% Inicialização do balanço harmônico

addpath('Tiro múltiplo');
tiro=tiro_multiplo;
tiro.din_func=@duffing;
tiro.din_jacobiana_x=@jacobiana_duffing_x;
tiro.din_jacobiana_pars=@jacobiana_duffing_pars;
        
tiro.autonomo=false;
tiro.complexo=false;
        
tiro.n_din=2;
tiro.p=4;

tiro.Nt=10;

tiro=tiro.montar_modelo();

%% 
% Numero de variaveis
modelo.n=tiro.n;
% Numero de parametros
modelo.p=tiro.p;
% Funcao residuo
modelo.residuo=@(u) residuo_duffing_tiro(tiro, u);
% Funcao jacobiana
modelo.jacobiana=@(u) jacobiana_duffing_tiro(tiro, u);
% Funcao de determinacao da estabilidade
modelo.estabilidade=@estabilidade_tiro_Floquet;
% Limites de interesse das variaveis             
limites_Q=repmat([-100, 100], tiro.n, 1);
        
limites_pars=[ -10,   10;
               -10,   10;
               -10,   10;
                 0,    3;];
modelo.limites=[limites_Q; limites_pars];
modelo.outdof1=modelo.n+modelo.p;
modelo.outdof2=3;
% pasta do modelo
modelo.folder='duffing_tiro';
modelo.data=tiro;

%% Parametros iniciais
x0=zeros(tiro.n_din, 1);
pars0=[0.01; 0.01; 1; 0.2];
X0=tiro.ponto_inicial_tiro(x0, pars0);

u0=[X0; pars0];
modelo.u0=u0;

H0=modelo.residuo(u0);
J0=modelo.jacobiana(u0);

h=1e-8;
for i=1:numel(u0)
    uf=u0;
    uf(i,1)=uf(i,1)+h;
    Hf=modelo.residuo(uf);
    Jf(:,i)=(Hf-H0)/h;
end

t0=vetor_tangente(J0);
estbilidade0=modelo.estabilidade(u0, J0, modelo);

%% parametros

% Passo e sentido iniciais
h0=0.0001;
w0=1;
       
% Numero de pontos maximo determinados
parametros.cont_max=10000;
% Numero maximo de iteracoes do corretor
parametros.i_max=6;
% Tamanho de passo minimo e maximo
parametros.h_min=1e-6;
parametros.h_max=1;
% Parametro de ajuste do passo
parametros.i_controle=3;
% Tolerancias
parametros.tol1=1e-6;
parametros.tol2=1e-6;

% Identificador do arquivo de saida de dados
id_saida=fopen(fullfile(modelo.folder, 'resultados_duffing_tiro1.txt'), 'w');

save(fullfile(modelo.folder, 'duffing_tiro.mat'), 'modelo', 'parametros');

end