classdef tiro_multiplo
    properties
        din_func %função do modelo dinâmico
        din_jacobiana_x %jacobiana do modelo dinâmico em relação às variáveis de estado
        din_jacobiana_pars %jacobiana do modelo dinâmico em relação aos parâmetros
        
        autonomo=false % verifica se o sistema dinâmico é autonomo (true ou false)
        complexo=false % para sistemas em que as variáveis de estado são pares complexo conjugados
        
        n_din % número de variáveis de estado do modelo dinâmico
        
        tol_abs=1e-8
        tol_rel=1e-8
        
        odefun=@ode45
        
        Nt %número de divisões do tempo
        
        n %número de variáveis método da continuação
        p %número de parâmetros do método da continuação
        
    end
    
    methods
        function obj=montar_modelo(obj)
            obj.n=obj.n_din*obj.Nt;
        end
        function R=residuo(obj, X, pars)
            global u_atual
            %Período
            T=1;
            
            if obj.autonomo==false
                fun=@(tt, xx) obj.din_func(tt, xx, pars);
            else
                %Para restrição de fase
                X_RF=u_atual(1:obj.n,:);
                pars_RF=u_atual(obj.n+1:obj.n+obj.p,:);
                R_RF=0;
                
                fun=@(tt, xx) obj.din_func(xx, pars);
                fun_RF=@(tt, xx) obj.din_func(xx, pars_RF);
            end
            
            options = odeset('RelTol',obj.tol_rel, 'AbsTol',obj.tol_abs);
            
            for i=1:obj.Nt
                
                x_i=X((i-1)*obj.n_din+1:obj.n_din*i,1);
                
                [~,Y]=ode45(fun, [T*(i-1)/obj.Nt T*i/obj.Nt], x_i, options);
                Y_i{i}=Y(end,:)';
                
            end
            for i=1:obj.Nt
                
                if i~=obj.Nt
                    intervalo1=[1:obj.n_din]+(obj.n_din*i);
                else
                    intervalo1=1:obj.n_din;
                end
                
                x_final=Y_i{i};
                
                if obj.autonomo==true
                    x_RF_i=X_RF(intervalo1,1);
                    F_RF_i=fun_RF(0, x_RF_i);
                    R_RF=R_RF+dot(x_RF_i-X(intervalo1), F_RF_i);
                end
                
                if(i~=obj.Nt)
                    X_final(i*obj.n_din+1:obj.n_din*(i+1),1)=x_final;
                else
                    X_final(1:obj.n_din,1)=x_final;
                end
            end
            %Cálculo do resíduo de lambuja
            
            R=X_final-X;
            if obj.autonomo==true
                R=[R;R_RF];
            end
        end
        
        function [R, Jx, Jp]=jacobiana(obj, X0, pars)
            global u_atual
            R=zeros(obj.n,1);
            R_RF=0;
            Jx=-eye(obj.n, obj.n);
            Jp=zeros(obj.n, obj.p);
            J_RF_x=zeros(1, obj.n);
            J_RF_p=zeros(1, obj.p);
            
            N_din=obj.n_din+obj.p; % número de variáveis totais do sistema expandido y=[x;pars]
            
            %Período
            T=1;
            
            if obj.autonomo==true
                %Para restrição de fase
                X_RF=u_atual(1:obj.n,:);
                pars_RF=u_atual(obj.n+1:obj.n+obj.p,:);
            
                fun_RF=@(tt, xx) obj.din_func(xx, pars_RF);
            end
            
            C0=eye(N_din, N_din);
            
            fun=@(tt, xx) din_monodromia(obj, tt, xx);
            
            options = odeset('RelTol',obj.tol_rel, 'AbsTol',obj.tol_abs);
            yC_i=cell(1,obj.Nt);
            
            for i=1:obj.Nt
                
                intervalo0=[1:obj.n_din]+(i-1)*obj.n_din;
                x_i=X0(intervalo0,1);
                y=[x_i;pars];
                
                
                yC0=[y;reshape(C0,[],1)];
                
                [~,yC]=ode45(fun, [T*(i-1)/obj.Nt T*i/obj.Nt], yC0, options);
                yC_i{i}=yC(end,:)';
                
            end
            
            for i=1:obj.Nt
                
                intervalo0=[1:obj.n_din]+(obj.n_din*(i-1));
                
                if i~=obj.Nt
                    intervalo1=[1:obj.n_din]+(obj.n_din*i);
                else
                    intervalo1=1:obj.n_din;
                end
                
                yC_final=yC_i{i};
                x_final_i=yC_final(1:obj.n_din,1);
                pars_final_i=yC_final(obj.n_din+1:obj.n_din+obj.p,1);
                
                
                J_i=reshape(yC_final(N_din+1:N_din+N_din*N_din, 1), N_din, N_din);
                Jx_i=J_i(1:obj.n_din, 1:obj.n_din);
                Jp_i=J_i(1:obj.n_din, obj.n_din+1:N_din);
                
                if obj.autonomo==true
                    x_RF_i=X_RF(intervalo1,1);
                    F_RF_i=fun_RF(0, x_RF_i);

                    J_RF_x_i=-F_RF_i';
                    
                    R_RF=R_RF+dot(x_RF_i-X0(intervalo1), F_RF_i);
                    J_RF_x(1, intervalo1)=J_RF_x_i;
                end
                
                R(intervalo1,1)=x_final_i-X0(intervalo1,1);    
                Jx(intervalo1, intervalo0)=Jx_i;
                Jp(intervalo1, :)=Jp_i;
            end
            
            if obj.autonomo==true
                R=[R;R_RF];
                Jx=[Jx; J_RF_x];
                Jp=[Jp; J_RF_p];
            end
        end
        
        function [X0]=ponto_inicial_tiro(obj,  x0, pars)
            
            if obj.autonomo==false
                fun=@(tt, xx) obj.din_func(tt, xx, pars);
            else
                fun=@(tt, xx) obj.din_func(xx, pars);
            end
            
            options = odeset('RelTol',obj.tol_rel, 'AbsTol',obj.tol_abs);
            
            for i=1:20
                [t,y]=ode45(fun, [0 1], x0, options);
                x0=y(end,:)';
            end
            
            [t,y]=ode45(fun, linspace(0, 1, obj.Nt+1), x0, options);
            X0=y(1:(end-1),:)';
            X0=reshape(X0, [], 1);
        end
        
    end
    methods(Access=private)
        function dFdt=din_monodromia(obj, t, yC)
            
            N_din=obj.n_din+obj.p; % número de variáveis totais do sistema expandido
            y=yC(1:N_din, :); %variáveis do sistema expandido [X; pars]
            C=yC(N_din+1:N_din+N_din*N_din,:); %segundo sistema expandido [u, vec(C)]
            
            C=reshape(C,N_din,N_din,[]);
            
            n_tempo=size(yC,2);
            
            X=y(1:obj.n_din,:);
            pars=y(obj.n_din+1:N_din,:);
            
            if obj.autonomo==false
                f=obj.din_func(t, X, pars);
                Jx=obj.din_jacobiana_x(t, X, pars);
                Jp=obj.din_jacobiana_pars(t, X, pars);
            else
                f=obj.din_func(X, pars);
                Jx=obj.din_jacobiana_x(X, pars);
                Jp=obj.din_jacobiana_pars(X, pars);
            end
            F=cat(1, f, zeros(obj.p,n_tempo));
            J=cat(1, cat(2, Jx, Jp), zeros(obj.p, N_din, n_tempo));
            
            M=reshape(pagemtimes(J, C), N_din*N_din, [], 1);
            
            dFdt=reshape([F;M], [],1);
            
        end
    end
end

