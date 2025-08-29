function Jp=jacobiana_pars(obj, Q, pars)

x=obj.Q_to_x(Q);

tempo=obj.tempo(1:end-1);

if obj.complexo==true
    x_aux=x(1:2:obj.M-1,:)+1i*x(2:2:obj.n_din,:);
    x(1:2:obj.n_din-1,:)=x_aux;
    x(2:2:obj.n_din,:)=conj(x_aux);
end

if obj.autonomo==false
    J=obj.din_jacobiana_pars(tempo, x, pars);
else
    J=obj.din_jacobiana_pars(x, pars);
end

if obj.complexo==true
    d_Cmplx_d_Real=kron(eye(obj.n_din/2, obj.n_din/2), [1, 1i; 1, -1i]);
    Jaux=pagemtimes(J,d_Cmplx_d_Real);
    J(1:2:obj.M-1,:,:)=real(Jaux(1:2:obj.n_din-1,:,:));
    J(2:2:obj.M,:,:)=imag(Jaux(1:2:obj.n_din-1,:,:));
end

J2 = reshape(permute(J, [1 3 2]), [], size(J, 2));
Jp=-obj.Einv*J2;

if obj.autonomo==true
    J_RF_pars=jacobiana_RF_pars(obj, Q, pars);
    Jp=[Jp;
       J_RF_pars];
end
    

end