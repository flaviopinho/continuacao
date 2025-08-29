function J=jacobiana_Q(obj, Q, pars)

x=obj.Q_to_x(Q);

tempo=obj.tempo(1:end-1);

if obj.complexo==true
    x_aux=x(1:2:obj.M-1,:)+1i*x(2:2:obj.n_din,:);
    x(1:2:obj.n_din-1,:)=x_aux;
    x(2:2:obj.n_din,:)=conj(x_aux);
end

if obj.autonomo==false
    J=obj.din_jacobiana_x(tempo, x, pars);
else
    J=obj.din_jacobiana_x(x, pars);
end

if obj.complexo==true
    d_Cmplx_d_Real=kron(eye(obj.n_din/2, obj.n_din/2), [1, 1i; 1, -1i]);
    Jaux=pagemtimes(J,d_Cmplx_d_Real);
    J(1:2:obj.M-1,:,:)=real(Jaux(1:2:obj.n_din-1,:,:));
    J(2:2:obj.M,:,:)=imag(Jaux(1:2:obj.n_din-1,:,:));
end


Kcell = cellfun(@sparse,  num2cell(J,[1,2]), 'uni',0  );
J=blkdiag(Kcell{:});
J_Q=obj.Einv*J*obj.E;

Jy_Q=2*pi*obj.Nabla;
J=Jy_Q-J_Q;

if(obj.autonomo==true)
    J_RF_Q=jacobiana_RF_Q(obj, Q, pars);
    J=[J;
       J_RF_Q];
end

end