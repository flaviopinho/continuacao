function plotar_resultados

% Script para plotar resultados de arquivos texto no formato especificado

% Configuração inicial
clear; clc; close all

%% Entrada de dados
arquivo = 'resultados_predador_presa*.txt';

col_x = 5;
col_y = 1;
col_z = 2;
dimensao = 3;

Titulo='Diagrama de bifurcação Modelo Predador-Presa';

%%

tamanho_ponto_regular=5;
tamanho_ponto_singular=15;

% Obter a lista de arquivos com o padrão especificado
files = dir(fullfile(fileparts(mfilename('fullpath')), arquivo));
if isempty(files)
    disp('Nenhum arquivo encontrado com o padrão especificado.');
    return;
end

% Inicializar variáveis para armazenar os dados
all_types = {};
all_data = [];

% Ler e combinar os dados de todos os arquivos
for i = 1:length(files)
    % Construir o caminho completo do arquivo
    filepath = fullfile(files(i).folder, files(i).name);
    
    % Ler o arquivo
    data = readtable(filepath, 'Delimiter', ';', 'ReadVariableNames', false);
    
    % Combinar os dados
    all_types = [all_types; data.Var1]; % Coluna 1 (tipo do ponto)
    num_columns = size(data, 2); % Número de colunas do arquivo
    all_data = [all_data; table2array(data(:, 2:num_columns))]; % Combinar colunas 2 em diante
end

% Número total de colunas disponíveis
num_columns = size(all_data, 2);

% Exibir opções de colunas para o usuário
disp('Colunas disponíveis para plotagem:');
for col = 1:num_columns
    fprintf('Coluna %d\n', col);
end

% Verificar se as colunas escolhidas são válidas
if col_x < 1 || col_x > num_columns || col_y < 1 || col_y > num_columns || ...
   (dimensao == 3 && (col_z < 1 || col_z > num_columns))
    error('Uma ou mais colunas escolhidas estão fora do intervalo válido. O script será encerrado.');
end

% Identificar tipos de pontos
stable_idx = strcmp(all_types, 'PRE'); % Estável
unstable_idx = strcmp(all_types, 'PRI'); % Instável
bifurcation_idx = strcmp(all_types, 'SN') | strcmp(all_types, 'H') | ...
                  strcmp(all_types, 'BPS') | strcmp(all_types, 'BPC');

% Selecionar os dados escolhidos pelo usuário
x = all_data(:, col_x);
y = all_data(:, col_y);
if dimensao == 3
    z = all_data(:, col_z);
end

if dimensao == 2
    % Plotar gráfico 2D
    figure;
    hold on;
    scatter(x(stable_idx), y(stable_idx), tamanho_ponto_regular, 'b', 'filled', 'DisplayName', 'Estável');
    scatter(x(unstable_idx), y(unstable_idx), tamanho_ponto_regular, 'r', 'filled', 'DisplayName', 'Instável');
    scatter(x(bifurcation_idx), y(bifurcation_idx), tamanho_ponto_singular, 'k', 'filled', 'DisplayName', 'Bifurcação/Ramificação');
    for i = find(bifurcation_idx)'
        text(x(i), y(i), [' ' all_types{i}], 'FontSize', 9, 'Color', 'k');
    end
    legend;
    xlabel(['Coluna ', num2str(col_x)]);
    ylabel(['Coluna ', num2str(col_y)]);
    title(Titulo);
    grid on;
    hold off;
    
elseif dimensao == 3
    % Plotar gráfico 3D
    figure;
    hold on;
    scatter3(x(stable_idx), y(stable_idx), z(stable_idx), tamanho_ponto_regular, 'b', 'filled', 'DisplayName', 'Estável');
    scatter3(x(unstable_idx), y(unstable_idx), z(unstable_idx), tamanho_ponto_regular, 'r', 'filled', 'DisplayName', 'Instável');
    scatter3(x(bifurcation_idx), y(bifurcation_idx), z(bifurcation_idx), tamanho_ponto_singular, 'k', 'filled', 'DisplayName', 'Bifurcação/Ramificação');
    % Exibir bifurcações com seus tipos
    for i = find(bifurcation_idx)'
        text(x(i), y(i), z(i), [' ' all_types{i}], 'FontSize', 9, 'Color', 'k');
    end
    legend;
    xlabel(['Coluna ', num2str(col_x)]);
    ylabel(['Coluna ', num2str(col_y)]);
    zlabel(['Coluna ', num2str(col_z)]);
    title(Titulo);
    grid on;
    view(3);
    hold off;
    
else
    disp('Dimensão inválida. O script será encerrado.');
end
