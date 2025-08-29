function Jp = jacobiana_duffing_pars(t, x, pars)
    % Extração dos parâmetros
    zeta = pars(1);    % Amortecimento normalizado
    alpha = pars(2);   % Rigidez cúbica
    F = pars(3);       % Amplitude do carregamento
    Omega = pars(4);   % Frequência da excitação
    
    % Fator de normalização no tempo
    T = 2 * pi / Omega;
    dT_dOmega = -2 * pi / Omega^2;  % Derivada de T em relação a Omega
    
    % Inicialização da matriz Jacobiana
    Jp = zeros(2, 4, length(t)); % Dimensão 3D para múltiplos valores de tempo
    
    % Primeira linha da Jacobiana (d(x1')/dp)
    Jp(1,1,:) = 0;
    Jp(1,2,:) = 0;
    Jp(1,3,:) = 0;
    Jp(1,4,:) = dT_dOmega * x(2,:); % Derivada de x1' em relação a Omega

    % Segunda linha da Jacobiana (d(x2')/dp)
    Jp(2,1,:) = T * (-2 * x(2,:));              % Derivada em relação a zeta
    Jp(2,2,:) = T * (-x(1,:).^3);               % Derivada em relação a alpha
    Jp(2,3,:) = T * cos(2 * pi * t(:).');       % Derivada em relação a F
    Jp(2,4,:) = dT_dOmega * (-2 * zeta * x(2,:) - x(1,:) - alpha * x(1,:).^3 + F * cos(2 * pi * t(:).')); % Derivada em relação a Omega
end