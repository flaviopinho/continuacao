%jacobiana_predador_presa.m
function J=jacobiana_predador_presa(u)

n=2;
np=3;

x=u(1:n,1);
p=u(n+1:n+np,1);

dFdx=[p(1)*(1 - x(1)) - p(1)*x(1) - x(2) - p(3)*p(2)*exp(-p(2)*x(1)), -x(1);
    p(1)*x(2), p(1)*x(1) - 1];
dFdp=[x(1)*(1 - x(1)), -p(3)*x(1)*exp(-p(2)*x(1)), -(1 - exp(-p(2)*x(1)));
    x(1)*x(2), 0, 0];

dgdx=[0, 0;
    0, 0];

dgdp=[1, 0, 0;
    0, 1, 0];

J=[dFdx, dFdp;
    dgdx, dgdp];

end