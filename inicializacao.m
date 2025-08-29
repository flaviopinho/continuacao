% inicializacao.m
function [u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao()

clc;
close('all');
fclose('all');

%addpath('predador_presa');
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_predador_presa();
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_ponto_de_ramificacao_predador_presa();

addpath('duffing');
[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_duffing();

%addpath('kerschen4');
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_kerschen4();
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_kerschen2();

%addpath('duffing_tiro');
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_duffing_tiro();

%addpath('kerschen4_tiro');
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_kerschen4_tiro();
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_kerschen2_tiro();