%ajustar_passo.m
function h=ajustar_passo(i, h, parametros)

h=max(min(h*min([max([sqrt(parametros.i_controle/i);0.5]);1.3]), parametros.h_max), parametros.h_min);

end
