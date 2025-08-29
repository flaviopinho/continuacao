function [estabilidade, tipo]=estabilidade_tiro_Floquet(u, J, modelo)

Nt=modelo.data.Nt;
n=modelo.data.n;
n_din=modelo.data.n_din;

for i=1:Nt
    if i~=Nt
        G{i}=J(i*n_din+1:n_din*(i+1),(i-1)*n_din+1:n_din*i);
    else
        G{i}=J(1:n_din,(i-1)*n_din+1:n_din*i);
    end
end
M=G{Nt};
for i=Nt-1:-1:1
    M=M*G{i};
end

mu=eig(M);

[r,ind]=max(abs(mu));

if(modelo.data.autonomo==false)
    estabilidade=r-1;
else
    estabilidade=r-(1+1e-4);
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