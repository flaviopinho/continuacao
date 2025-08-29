function J=jacobiana_duffing_tiro(tiro, u)

n=tiro.n; % número de variáveis do método do balanço harmônico
p=tiro.p; % número de parâmetros

X=u(1:n,1); 
pars=u(n+1:n+p,1);

[~,dF_dX,dF_dpars]=jacobiana(tiro, X, pars);

dg_dX=zeros(3, n);
dg_dpars=[1, 0, 0, 0;
       0, 1, 0, 0;
       0, 0, 1, 0]; 

J=[dF_dX, dF_dpars;
   dg_dX, dg_dpars];

end