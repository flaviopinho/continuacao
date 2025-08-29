function y = duffing(t, x, pars)
    % Extração dos parâmetros
    zeta = pars(1);    % Amortecimento normalizado
    alpha = pars(2);   % Rigidez cúbica
    F = pars(3);       % Amplitude do carregamento
    Omega = pars(4);   % Frequência da excitação
    
    % Inicialização da matriz de saída
    y = zeros(size(x));  % Mesmas dimensões de x

    % Cálculo para cada estado (coluna de x)
    T=2*pi/Omega;
    y(1, :) = T*x(2, :);
    y(2, :) = T*(-2 * zeta * x(2, :) - x(1, :) - alpha * x(1, :).^3 + F * cos(2* pi * t(:).'));
end
