%vetor_tangente.m
function [t,f]=vetor_tangente(J)

[n,m]=size(J);

[Q,R]=qr(J');
f=det(Q)*det(R(1:n,:));
t=sign(f)*Q(:, m);

end