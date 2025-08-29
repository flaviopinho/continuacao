function J=jacobiana_kerschen4_bh_autonomo(bh, u)

n=bh.n; % número de variáveis do método do balanço harmônico
p=bh.p; % número de parâmetros

Q=u(1:n,1); 
pars=u(n+1:n+p,1);

dF_dQ=jacobiana_Q(bh, Q, pars);
dF_dpars=jacobiana_pars(bh, Q, pars);
dg_dQ=zeros(1, n);
dg_dpars=[1, 0, 0]; 

J=[dF_dQ, dF_dpars;
   dg_dQ, dg_dpars];

end