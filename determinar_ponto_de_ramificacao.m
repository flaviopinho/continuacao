% determinar_ponto_de_ramificacao.m
function determinar_ponto_de_ramificacao(u, v, modelo, parametros, id_saida)

u1=u;
u2=v;

[J1]=feval(modelo.jacobiana,u1);
[J2]=feval(modelo.jacobiana,u2);

[~, f_1]=vetor_tangente(J1);
[~, f_2]=vetor_tangente(J2);


for j=1:10
    um=u1-f_1/(f_2-f_1)*(u2-u1);

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
    
    [~, f_m]=vetor_tangente(J_um);
    
    if(abs(f_m)<parametros.tol1*0.01)
        break;
    end
    
    if(f_1*f_m<0)
        u2=um;
        f_2=f_m;
    elseif(f_2*f_m<0)
        u1=um;
        f_1=f_m;
    end
    
end

if(size(J_um,2)-rank(J_um, parametros.tol1)==2)
    tipo='BPS'; % Ponto de ramificacao simples
else
    tipo='BPC'; % Ponto de ramificacao complexo
end

fprintf(id_saida, '%s;', tipo);
for i=1:numel(um)
    fprintf(id_saida, ' %e;', um(i,1));
end
fprintf(id_saida, '\n');
hold on;
plot(um(end,1), um(1,1), 'ko');

end