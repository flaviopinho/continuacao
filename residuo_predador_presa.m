%residuo_predador_presa.m
function H=residuo_predador_presa(u)

n=2;
np=3;

x=u(1:n,1);
p=u(n+1:n+np,1);

F=[p(1)*x(1)*(1-x(1))-x(1)*x(2)-p(3)*(1-exp(-p(2)*x(1)));
    -x(2)+p(1)*x(1)*x(2)];

g=[p(1)-3
    p(2)-5];

H=[F;
    g];

end