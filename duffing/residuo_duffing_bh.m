function R=residuo_duffing_bh(bh, u)

n=bh.n; % número de variáveis do método do balanço harmônico
p=bh.p; % número de parâmetros

Q=u(1:n,1); 
pars=u(n+1:n+p,1);

% Extração dos parâmetros
zeta = pars(1);    % Amortecimento normalizado
alpha = pars(2);   % Rigidez cúbica
F = pars(3);       % Amplitude do carregamento
Omega = pars(4);   % Frequência da excitação

Fr=residuo_Q(bh, Q, pars);
g=[zeta-pars(1); alpha-pars(2); F-pars(3)];
R=[Fr; g];

end