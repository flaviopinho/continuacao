function Mon=monodromia(obj, u)
Q=u(1:obj.n, :);
pars=u(obj.n+1:(obj.n+obj.p), :);
x=Q_to_x(obj, Q);
x0=x(:,1);
M=eye(obj.n_din, obj.n_din);

y0=[x0;reshape(M, [], 1)];
func2=@(tt, xx) FunMon(obj, tt, xx, pars);
[tt,xx]=ode45(func2, [0 1], y0);
xx=xx';
y1=xx(:,end);
x1=xx(1:obj.n_din,1);
Mon=reshape(y1(obj.n_din+1:end,1), obj.n_din, obj.n_din);
end

function F=FunMon(obj, t, y, pars)
x=y(1:obj.n_din);

if obj.complexo==true
    x_aux=x(1:2:obj.M-1,:)+1i*x(2:2:obj.n_din,:);
    x(1:2:obj.n_din-1,:)=x_aux;
    x(2:2:obj.n_din,:)=conj(x_aux);
end

M=reshape(y(obj.n_din+1:end,1), obj.n_din, obj.n_din);
if obj.autonomo==true
    F=obj.din_func(x, pars);
    J=obj.din_jacobiana_x(x, pars)*M;
else
    F=obj.din_func(t, x, pars);
    J=obj.din_jacobiana_x(t, x, pars)*M;
end

if obj.complexo==true
    F_real(1:2:obj.M-1,:)=real(F(1:2:obj.M-1,:));
    F_real(2:2:obj.M,:)=imag(F(1:2:obj.M-1,:));
    F=F_real;
    
    d_Cmplx_d_Real=kron(eye(obj.n_din/2, obj.n_din/2), [1, 1i; 1, -1i]);
    Jaux=pagemtimes(J,d_Cmplx_d_Real);
    J(1:2:obj.M-1,:,:)=real(Jaux(1:2:obj.n_din-1,:,:));
    J(2:2:obj.M,:,:)=imag(Jaux(1:2:obj.n_din-1,:,:)); 
end

F=[F;reshape(J, [], 1)];
end