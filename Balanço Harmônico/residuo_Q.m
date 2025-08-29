function R=residuo_Q(obj, Q, pars)

x=obj.Q_to_x(Q);

tempo=obj.tempo(1:end-1);

if (obj.complexo==true)
    x_real=x;
    x_aux=x_real(1:2:obj.M-1,:)+1i*x_real(2:2:obj.M,:);
    x_cmplx(1:2:obj.M-1,:)=x_aux;
    x_cmplx(2:2:obj.M,:)=conj(x_aux);
    x=x_cmplx;
end

if obj.autonomo==false
    F=obj.din_func(tempo, x, pars);
else
    F=obj.din_func(x, pars);
end

if (obj.complexo==true)
    F_real(1:2:obj.M-1,:)=real(F(1:2:obj.M-1,:));
    F_real(2:2:obj.M,:)=imag(F(1:2:obj.M-1,:));
    F=F_real;
end

F_Q=obj.Einv*reshape(F, obj.n_din*obj.Nt, 1, 1);
R=2*pi*obj.Nabla*Q-F_Q;

if obj.autonomo==true
    RF=residuo_RF(obj, Q, pars);
    R=[R;
       RF];
end