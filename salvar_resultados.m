function salvar_resultados(u, teste, id_saida, modelo)

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
        plot(u(modelo.outdof1,1), u(modelo.outdof2,1), 'b.');
    else
        plot(u(modelo.outdof1,1), u(modelo.outdof2,1), 'r.');
    end
    drawnow;
end