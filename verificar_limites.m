%verificar_limites.m
function teste=verificar_limites(u, modelo)
    teste=false;
    
    u_min=modelo.limites(:,1);
    u_max=modelo.limites(:,2);
    
    if(min(u_max-u)<0)
        teste=true;
    elseif(min(u-u_min)<0)
        teste=true;
    end
end