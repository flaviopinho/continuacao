function dxdt = kerschen4(t, x, pars)
    % Extração dos parâmetros
    F = pars(1);       % Amplitude do carregamento
    c = pars(2);       % Amortecimento
    Omega = pars(3);   % Frequência da excitação
    
    % Inicialização da matriz de saída
    dxdt = zeros(size(x));  % Mesmas dimensões de x

    % Cálculo para cada estado (coluna de x)
    T=2*pi/Omega;
    % Equações do sistema
    dxdt(1,:) = T*x(3,:);
    dxdt(2,:) = T*x(4,:);
    dxdt(3,:) = T*(-2*c*x(3,:) + c*x(4,:) - 2*x(1,:) + x(2,:) - 0.5*x(1,:).^3 + F*cos(2*pi*t(:).'));
    dxdt(4,:) = T*(-2*c*x(4,:) + c*x(3,:) - 2*x(2,:) + x(1,:));
end