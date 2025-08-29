function Jx = jacobiana_duffing_x(t, x, pars)
    % Extração dos parâmetros
    zeta = pars(1);    % Amortecimento normalizado
    alpha = pars(2);   % Rigidez cúbica
    F = pars(3);       % Amplitude do carregamento
    Omega = pars(4);   % Frequência da excitação

    % Fator de normalização no tempo
    T = 2 * pi / Omega;

    % Inicialização da matriz Jacobiana
    Jx = zeros(2, 2, size(x,2)); % Dimensão 3D para múltiplos estados

    % Primeira linha da Jacobiana (d(x1')/dx)
    Jx(1,1,:) = 0;
    Jx(1,2,:) = T;

    % Segunda linha da Jacobiana (d(x2')/dx)
    Jx(2,1,:) = T * (-1 - 3 * alpha * x(1,:).^2); % Derivada em relação a x1
    Jx(2,2,:) = T * (-2 * zeta);                 % Derivada em relação a x2
end