function J = jacobiana_kerschen4_x(t, x, pars)
    % Extração dos parâmetros
    F = pars(1);       % Amplitude do carregamento
    c = pars(2);       % Amortecimento
    Omega = pars(3);   % Frequência da excitação
    T = 2 * pi / Omega;
    
    % Inicialização da matriz Jacobiana
    J = zeros(4, 4, size(x, 2));
    
    % Derivadas parciais
    J(1,3,:) = T;
    J(2,4,:) = T;
    J(3,1,:) = T * (-2 - 1.5*x(1,:).^2);
    J(3,2,:) = T;
    J(3,3,:) = T * (-2*c);
    J(3,4,:) = T * (c);
    J(4,1,:) = T;
    J(4,2,:) = T * (-2);
    J(4,3,:) = T * (c);
    J(4,4,:) = T * (-2*c);
end