%continuacao.m
function continuacao()

%Inicializacao das variaveis iniciais, do modelo e dos parametros
[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao();
%[u, t_u, estabilidade_u, h, w, modelo, parametros, id_saida]=inicializacao_ponto_de_ramificacao(); 

% Algoritmo de continuacao
for cont=1:parametros.cont_max
    
    if(h<parametros.h_min)
        fprintf('-------- Tamanho do passo menor que o minimo. --------\n');
        break;
    end
    
    % preditor
    v=u+h*w*t_u;

    for i=1:parametros.i_max
        %Corretor
        [H]=feval(modelo.residuo,v);
        [J_v]=feval(modelo.jacobiana,v);
        dv=pinv(J_v)*H;
        v=v-dv;

        erro1=norm(H, 2);
        erro2=norm(dv, 2)/norm(v,2);
        
        if(erro1<parametros.tol1 || erro2<parametros.tol2)
            break;
        end

        if(i==parametros.i_max)
            % reducao de passo caso numero de iteracoes do corretor ultrapasse o limite
            v=u;
            h=h/2;
            fprintf('-------- Corretor nao convergiu. Reducao de passo. --------\n');
            continue;
        end
    end

    % Verificar bifurcacoes
    [estabilidade_v]=feval(modelo.estabilidade, v, J_v, modelo);
    
    if(estabilidade_u*estabilidade_v<0)
        fprintf('-------- Bifurcacao. --------\n');
        determinar_ponto_de_bifurcacao(u, v, modelo, parametros, id_saida);
    end
    
    % Verificar ponto de ramificacao
    t_v=vetor_tangente(J_v);
    
    if(dot(t_u,t_v)<0)
        fprintf('-------- Ponto de ramificacao. --------\n');
        determinar_ponto_de_ramificacao(u, v, modelo, parametros, id_saida);
        w=-w;
    end

    %Inicializacao das variaveis para proximo ponto da curva
    u=v;
    t_u=t_v;
    J_u=J_v;
    estabilidade_u=estabilidade_v;
    
    % Impressão de resultados parciais em tela
    fprintf('n=%i   i=%i  erro1=%e   erro2=%e   h=%e  estabilidade=%e \n', cont, i, erro1, erro2, h, estabilidade_v);
    
    % Saida dos resultados
    salvar_resultados(u, estabilidade_u, id_saida);
    
    %Ajustar tamanho do passo
    [h]=ajustar_passo(i, h, parametros);
    
    % Verificar intevalos das variaveis para finalizacao do loop
    if(verificar_limites(u, modelo))
        break;
    end
end

fclose('all');

end