function RF=residuo_RF(obj, Q, pars)
    RF=2*pi*obj.E(obj.RF_var,:)*obj.Nabla*Q;
end