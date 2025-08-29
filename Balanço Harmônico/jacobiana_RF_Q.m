function J_RF=jacobiana_RF_Q(obj, Q, pars)
    J_RF=2*pi*obj.E(obj.RF_var,:)*obj.Nabla;

end