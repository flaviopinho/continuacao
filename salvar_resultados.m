function salvar_resultados(u, teste, id_saida)

    if(teste<0)
        tipo='PRE'; % Ponto regular estável
    else
        tipo='PRI'; % Ponto regular instável
    end
    fprintf(id_saida, '%s;', tipo);
    for i=1:numel(u)
        fprintf(id_saida, ' %e;', u(i,1));
    end
    fprintf(id_saida, '\n');
    hold on;
    if(teste<0)
        plot(u(end,1), u(1,1), 'b.');
    else
        plot(u(end,1), u(1,1), 'r.');
    end
end