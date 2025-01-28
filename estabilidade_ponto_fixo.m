%estabilidade_ponto_fixo.m
function [estabilidade, tipo]=estabilidade_ponto_fixo(u, J, modelo)

Jx=J(1:modelo.n, 1:modelo.n);

% Determinacao dos autovalores de Jx
[~, auto_valores] = eig(Jx);

% Partes reais e imaginarias dos autovalores
parte_real = real(diag(auto_valores));
parte_imaginaria = imag(diag(auto_valores));

estabilidade=max(parte_real);

if nargout>1
    % Analise do tipo de bifurcacao
    index_real_positivo=parte_real > 0;
    index_real_negativo=parte_real < 0;
    num_real_positivo = sum(index_real_positivo);
    num_real_negativo = sum(index_real_negativo);
    
    if num_real_positivo == 2 && any(parte_imaginaria(index_real_positivo) ~= 0)
        tipo='H'; % Hopf 
    elseif num_real_positivo == 1 && num_real_negativo >= 0
        tipo='SN'; % Ponto de sela
    elseif num_real_positivo == 0 && num_real_negativo >= 0
        tipo='PR'; % Ponto regular
    else
        tipo='BC'; % Bifurcacao complexa
    end
end

end