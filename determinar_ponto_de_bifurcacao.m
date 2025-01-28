% determinar_ponto_de_bifurcacao.m
function determinar_ponto_de_bifurcacao(u, v, modelo, parametros, id_saida)

u1=u;
u2=v;

J1=feval(modelo.jacobiana,u1);
J2=feval(modelo.jacobiana,u2);

[estabilidade_u1, tipo1]=feval(modelo.estabilidade, u1, J1, modelo);
[estabilidade_u2, tipo2]=feval(modelo.estabilidade, u2, J2, modelo);

if(estabilidade_u1>0)
    tipo=tipo1;
else
    tipo=tipo2;
end

for j=1:10
    um=u1-estabilidade_u1/(estabilidade_u2-estabilidade_u1)*(u2-u1);

    for i=1:parametros.i_max
        %Corretor
        [H]=feval(modelo.residuo,um);
        [J_um]=feval(modelo.jacobiana,um);
        dum=lsqminnorm(J_um,H); % dum=pinv(J_um)*H;
        um=um-dum;
        
        erro1=norm(H, 2);
        erro2=norm(dum, 2)/norm(um,2);
        
        if(erro1<parametros.tol1 || erro2<parametros.tol2)
            break;
        end
    end
    
    [estabilidade_um]=feval(modelo.estabilidade, um, J_um, modelo);
    
    if(abs(estabilidade_um)<parametros.tol1)
        break;
    end
    
    if(estabilidade_u1*estabilidade_um<0)
        u2=um;
        estabilidade_u2=estabilidade_um;
    elseif(estabilidade_u2*estabilidade_um<0)
        u1=um;
        estabilidade_u1=estabilidade_um;
    end
    
end

fprintf(id_saida, '%s;', tipo);
for i=1:numel(um)
    fprintf(id_saida, ' %e;', um(i,1));
end
fprintf(id_saida, '\n');
hold on;
plot(um(end,1), um(1,1), 'ks');

end