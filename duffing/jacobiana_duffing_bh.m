function J=jacobiana_duffing_bh(bh, u)

n=bh.n; % número de variáveis do método do balanço harmônico
p=bh.p; % número de parâmetros

Q=u(1:n,1); 
pars=u(n+1:n+p,1);

% Extração dos parâmetros
zeta = pars(1);    % Amortecimento normalizado
alpha = pars(2);   % Rigidez cúbica
F = pars(3);       % Amplitude do carregamento
Omega = pars(4);   % Frequência da excitação

dF_dQ=jacobiana_Q(bh, Q, pars);
dF_dpars=jacobiana_pars(bh, Q, pars);
dg_dQ=zeros(3, n);
dg_dpars=[1, 0, 0, 0;
       0, 1, 0, 0;
       0, 0, 1, 0]; 

J=[dF_dQ, dF_dpars;
   dg_dQ, dg_dpars];

end