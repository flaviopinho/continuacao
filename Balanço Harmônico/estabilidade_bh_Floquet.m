function [estabilidade, tipo]=estabilidade_bh_Floquet(u, J, modelo)

M=monodromia(modelo.data, u);
mu=eig(M);

[r,ind]=max(abs(mu));

if(modelo.data.autonomo==false)
    estabilidade=r-1;
else
    estabilidade=r-(1+1e-6);
end

if(nargout>1)
    mu=mu(ind);
    if estabilidade<0
        tipo="PR";
    else
        % Identify and store type of BP
        if abs(angle(mu))<1e-3
            % Critical Floquet multiplier appears to leave unit circle 
            % via +1. This indicates a turning point or branching
            % bifuraction. In this case example it is a turning point (TP).
            tipo="SN";
        elseif abs(abs(angle(mu))-pi)<1e-3
            % Critical Floquet multiplier appears to leave unit circle 
            % via -1. This indicates a period doubling bifurcation (PD).
            tipo="PD";
        else
            % Critical Floquet multiplier appears to leave unit circle as
            % complex-valued pair. This indicates a Neimark-Sacker
            % bifurcation (NS).
            tipo="NS";
        end
    end
end