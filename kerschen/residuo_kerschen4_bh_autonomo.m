function R=residuo_kerschen4_bh_autonomo(bh, u)

n=bh.n; % número de variáveis do método do balanço harmônico
p=bh.p; % número de parâmetros

Q=u(1:n,1); 
pars=u(n+1:n+p,1);

% Extração dos parâmetros
F = pars(1);       % Amplitude do carregamento
c = pars(2);       % Amortecimento
Omega = pars(3);   % Frequência da excitação

Fr=residuo_Q(bh, Q, pars);
g=[F-pars(1)];
R=[Fr; g];

end