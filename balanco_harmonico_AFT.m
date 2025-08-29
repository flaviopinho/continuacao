% 
%% Variáveis importantes utilizadas no código
% *u*=[ *Q* , *pars* ] : variáveis do método da continuação
%
% *Q* : vetor com coeficientes de Fourier
%
% *pars* : vetor de parâmetros
%
% *x* : vetor das variáveis de estado
%
% * no caso complexo, esse vetor contem as partes reais e imginárias das
% variáveis de estado
%
% * no caso real, esse vetor contem os deslocamentos e velocidades
% 
% *X* : vetor das variáveis de estado para todos os instantes de tempo:
%
% |Nt| : número de instantes de tempo
%
% |H| : número de harmônicos
%
% |tempo=[0, 1/Nt, ..., (Nt-1)/Nt, 1]|


%% Definição da função din_func
% Para sistemas não autônomos:
%
% $$\mathbf{y}=\mathbf{f}(\tau, \mathbf{x}, \mathbf{pars})$$
%
% Para sistemas autônomos:
%
% $$\mathbf{y}'=\mathbf{f}(\mathbf{x}, \mathbf{pars})$$
%
% tal que:
%
% $$ 0\leq\tau=t/T \leq 1$$
%
% $$ T=2\pi/\Omega $$
%
% OBs.: A função din_func DEVE ser vetorizável
%
%% Definição da funções din_jacobiana_x e din_jacobiana_pars
%
% $$\mathbf{Jx}=\partial\mathbf{f}/\partial \mathbf{x} $$
%
% $$\mathbf{Jp}=\partial\mathbf{f}/\partial \mathbf{pars} $$
%
% OBs.: As funções din_jacobiana_x e din_jacobiana_pars DEVEM ser
% vetorizáveis
%
%% Função resíduo pelo método do balanço harmônico
%
% $$ \mathbf{y}-\mathbf{f}(\tau,\mathbf{x},\mathbf{pars})=\mathbf{0} $$
%
% $$ 2\pi \nabla \cdot  \mathbf{Q}(\mathbf{X})-\mathbf{Einv} \cdot \mathbf{f}(\mathrm{tempo}, \mathbf{X}(\mathbf{Q}), \mathbf{pars})=\mathbf{0} $$
%% Restrição de fase
%
% $$ \mathrm{RF}=\mathbf{y}[1](0)=\mathbf{0}$$
%
% $$ \mathrm{RF}=2\pi\mathbf{E}[1,:] \cdot \nabla \cdot \mathbf{Q} = \mathbf{0} $$
%
classdef balanco_harmonico_AFT
    properties(Access=public)
        
        din_func %função do modelo dinâmico
        din_jacobiana_x %jacobiana do modelo dinâmico em relação às variáveis de estado
        din_jacobiana_pars %jacobiana do modelo dinâmico em relação aos parâmetros
        
        autonomo=false % verifica se o sistema dinâmico é autonomo (true ou false)
        complexo=false % para sistemas em que as variáveis de estado são pares complexo conjugados
        
        n_din % número de variáveis de estado do modelo dinâmico
        
        H %número de harmônicos
        Nt %número de divisões do tempo
        
        n %número de variáveis método da continuação
        p %número de parâmetros do método da continuação
        
        tempo % instantes de tempo discretos
    
        E % Operador utilizado para sair de Q  para U, ou seja, do domínio da frequÊncia para o domínio do tempo
        Einv % Operador utilizado para sair de U  para Q, ou seja, do domínio do tempo para o domínio da frequência
        
        nabla % Matriz utilizada para o cálculo de Nabla
        Nabla % Operador utilizado para cálculo das velocidades
        
        RF_var=1 % restrição de fase
        
    end
    
    methods(Access=public)
        
        %% Relação entre Q e x
        %
        % $$ \mathbf{X}=[\mathbf{x}(0), \mathbf{x}(1/N_t), \dots,
        % \mathbf{x}((N_t-1)/N_t)] $$
        %
        % $$ \mathbf{Q}=\mathbf{Einv} \cdot \mathrm{vec}(\mathbf{y}) $$
        %
        % $$ \mathrm{vec}(\mathbf{X})=\mathbf{E} \cdot \mathbf{Q} $$
        %
        function x=Q_to_x(obj, Q)
            x=obj.E*Q;
            x=reshape(x, obj.n_din, obj.Nt);
        end
        
        function Q=x_to_Q(obj, x)
            x_vec=reshape(x, obj.n_din*obj.Nt, 1);
            Q=obj.Einv*x_vec;
        end
        %% Relação entre y e Q
        %
        % $$ \mathbf{Y}=[\mathbf{y}(0), \mathbf{y}(1/N_t), \dots,
        % \mathbf{y}((N_t-1)/N_t)] $$
        %
        % $$ \mathrm{vec}(\mathbf{Y})=2\pi\mathbf{E} \cdot \nabla \cdot \mathbf{Q} $$
        %
        function v=Q_to_v(obj, Q)
            v=obj.E*2*pi*obj.Nabla*Q;
            v=reshape(v, obj.n_din, obj.N);
        end
        
        function obj=montar_modelo(obj)
            obj.tempo=linspace(0, 1, obj.Nt+1);
            obj.E=obj.MatrizE();
            obj.Einv=obj.MatrizEInv();
            obj.nabla=obj.Matriz_nabla();
            obj.Nabla=kron(obj.nabla, speye(obj.n_din, obj.n_din));
            
            obj.n=obj.n_din*(2*obj.H+1);
        end
        
    end
    
    methods(Access=private)
        function E=MatrizE(obj)
            E1=zeros(obj.Nt, (2*obj.H+1));
            E1(:,1)=ones(obj.Nt, 1);
            for i=1:obj.H
                for j=1:obj.Nt
                    E1(j, 2*i)=cos(2*pi*(j-1)*i/obj.Nt);
                    E1(j, 2*i+1)=sin(2*pi*(j-1)*i/obj.Nt);
                end
            end
            E=kron(E1, speye(obj.n_din, obj.n_din));
        end
        
        function Einv=MatrizEInv(obj)
            Einv1=zeros((2*obj.H+1), obj.Nt);
            Einv1(1,:)=1/2*ones(obj.Nt, 1);
            for i=1:obj.H
                for j=1:obj.Nt
                    Einv1(2*i, j)  =cos(2*pi*i*(j-1)/obj.Nt);
                    Einv1(2*i+1, j)=sin(2*pi*i*(j-1)/obj.Nt);
                end
            end
            Einv=kron(2/obj.Nt*Einv1, speye(obj.n_din, obj.n_din));
        end
        
        function nabla=Matriz_nabla(obj)
            NablaK{1}=0;
            for k=1:obj.H
                NablaK{k+1}=[0 k; -k 0];
            end
            nabla=sparse(blkdiag(NablaK{:}));
        end
        
    end
    
end