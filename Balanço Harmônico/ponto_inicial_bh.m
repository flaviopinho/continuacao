function Q0=ponto_inicial_bh(bh, x0, pars0)

if bh.autonomo==false
    fun=@(tt,xx) bh.din_func(tt, xx, pars0);
else
    fun=@(tt,xx) bh.din_func(xx, pars0);
end

for i=1:20
    [t,y]=ode45(fun, bh.tempo, x0);
    x0=y(end,:)';
end

x0=y(1:end-1, :)';
Q0=bh.x_to_Q(x0);

end