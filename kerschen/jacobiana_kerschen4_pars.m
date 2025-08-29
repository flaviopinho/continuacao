function Jp = jacobiana_kerschen4_pars(t, x, pars)
    % Extração dos parâmetros
    F = pars(1);       % Amplitude do carregamento
    c = pars(2);   % Frequência da excitação
    Omega = pars(3);   % Frequência da excitação
    T = 2 * pi / Omega;
    
    % Inicialização da matriz Jacobiana em relação aos parâmetros
    Jp = zeros(4, 2, size(x, 2));
    
    % Derivadas parciais em relação aos parâmetros
    Jp(3,1,:) = T * cos(2 * pi * t(:).'); % Parcial de dx3/dF
    Jp(3,2,:) = T * (-2*x(3,:) + x(4,:)); % Parcial de dx3/dc
    Jp(4,2,:) = T * (-2*x(4,:) + x(3,:)); % Parcial de dx4/dc
    Jp(3,3,:) = -(2 * pi / Omega^2) * (-2*c*x(3,:) + c*x(4,:) - 2*x(1,:) + x(2,:) - 0.5*x(1,:).^3 + F*cos(2*pi*t(:).'));
    Jp(1,3,:) = -(2 * pi / Omega^2) * x(3,:); % Parcial de dx1/dOmega
    Jp(2,3,:) = -(2 * pi / Omega^2) * x(4,:); % Parcial de dx2/dOmega
    Jp(4,3,:) = -(2 * pi / Omega^2) * (-2*c*x(4,:) + c*x(3,:) - 2*x(2,:) + x(1,:)); % Parcial de dx4/dOmega
end