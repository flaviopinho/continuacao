%determinar_novas_ramificacoes.m
function [t1, t2]=determinar_novas_ramificacoes(J, tol=1e-6)

M=null_tol(J, tol);
if(size(M,2)~=2)
    error('Matriz J nao representa um ponto de ramificacao simples').
end
t1=M(:,1);
t2=M(:,2);

end