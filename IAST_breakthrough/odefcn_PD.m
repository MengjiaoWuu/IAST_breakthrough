function dydt=odefcn_PD(t,y,v,P,T,eps,rho,Dax,Nz,delta_z,isotherm_data_type,component_info,r_p,mu)
    % parameters
    Nc=size(component_info,1);
    y_if=vertcat(component_info.y_if);
    k_i=vertcat(component_info.k_i);
    
    interpolation_method='spline';
    R=8.3144598;% gas constant (m^3 Pa/(mol k))
    Nf=1+2*Nc;
    dydt=zeros(Nf*Nz,1);
    Delta=10^-(10);


    for j = 1 : Nz 
        if j==1

            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            P_j=y(Nf*(j-1)+2*Nc+1);
            
            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);
            P_jp1=y(Nf*j+2*Nc+1);


            %j+0.5
            y_i_jph=y_i_j;
            P_jph=P_j;
            v_jph=(-4/150)*((eps/(1-eps))^2)*(r_p^2/mu)*(1/delta_z)*(P_jp1-P_j);

            %j-0.5=1-0.5=0.5
            v_jmh=v;
            y_i_jmh=(Dax.*y_i_j+v_jmh.*(delta_z./2).*y_if)./(v_jmh.*(delta_z./2)+Dax);
            P_jmh=P_j+(delta_z/2)*(150/4)*(1/r_p^2)*((1-eps)/eps)^2*mu*v_jmh;
            % equilibrium uptake
            P_i_j=y_i_j*P_j;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(component_info,P_i_j,Nc).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(component_info,P_i_j).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=component_info(i).Fs(P_i_j').*rho;
                end
            end 

            % ode 
            sumdqdt=0;
            for i=1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                sumdqdt=sumdqdt+dydt(Nf*(j-1)+Nc+i);
            end

             dydt(Nf*(j-1)+2*Nc+1)=(-1/delta_z)*(P_jph*v_jph-P_jmh*v_jmh)-((1-eps)/eps)*R*T*sumdqdt;% total mass balance

            for i = 1: Nc
                dydt(Nf*(j-1)+i)=(1/P_j)*Dax*(1/delta_z)*(P_jph*(y_i_jp1(i)-y_i_j(i))/delta_z-P_jmh*(y_i_j(i)-y_i_jmh(i))/(delta_z/2))...
                -(1/P_j)*(1/delta_z)*(y_i_jph(i)*P_jph*v_jph-y_i_jmh(i)*P_jmh*v_jmh)...
                -(R*T/P_j)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i)-(y_i_j(i)/P_j)*dydt(Nf*(j-1)+2*Nc+1);% gas phase mass balance
            end 
                
        elseif j==2

            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            P_j=y(Nf*(j-1)+2*Nc+1);
            
            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);
            P_jm1=y(Nf*(j-2)+2*Nc+1);

            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);
            P_jp1=y(Nf*j+2*Nc+1);


            %j+0.5
            y_i_jph=y_i_j;
            P_jph=P_j;
            v_jph=(-4/150)*((eps/(1-eps))^2)*(r_p^2/mu)*(1/delta_z)*(P_jp1-P_j);

            %j-0.5
            v_jmh=(-4/150)*((eps/(1-eps))^2)*((r_p^2)/mu)*(1/delta_z)*(P_j-P_jm1);
            y_i_jmh=y_i_jm1;
            P_jmh=P_jm1;
            
            % equilibrium uptake
            P_i_j=y_i_j*P_j;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(component_info,P_i_j,Nc).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(component_info,P_i_j).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=component_info(i).Fs(P_i_j').*rho;
                end
            end 
            % ode 
            sumdqdt=0;
            for i=1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                sumdqdt=sumdqdt+dydt(Nf*(j-1)+Nc+i);
            end

             dydt(Nf*(j-1)+2*Nc+1)=(-1/delta_z)*(P_jph*v_jph-P_jmh*v_jmh)-((1-eps)/eps)*R*T*sumdqdt;% total mass balance

            for i = 1: Nc
                dydt(Nf*(j-1)+i)=(1/P_j)*Dax*(1/delta_z)*(P_jph*(y_i_jp1(i)-y_i_j(i))/delta_z-P_jmh*(y_i_j(i)-y_i_jm1(i))/delta_z)...
                -(1/P_j)*(1/delta_z)*(y_i_jph(i)*P_jph*v_jph-y_i_jmh(i)*P_jmh*v_jmh)...
                -(R*T/P_j)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i)-(y_i_j(i)/P_j)*dydt(Nf*(j-1)+2*Nc+1);% gas phase mass balance
            end 
           

        elseif j==Nz
            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            P_j=P;

            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);
            P_jm1=y(Nf*(j-2)+2*Nc+1);

            %j+0.5=Nz+0.5(end)
            y_i_jph=y_i_j;
            P_jph=P_j;
            v_jph=(-4/150)*(r_p^2)*((eps/(1-eps))^2)*((r_p^2)/mu)*(2/delta_z)*(P_jph-P_j);

            %j-0.5
            y_i_jmh=y_i_jm1;
            P_jmh=P_jm1;
            v_jmh=(-4/150)*((eps/(1-eps))^2)*((r_p^2)/mu)*(1/delta_z)*(P_j-P_jm1);
 
            % equilibrium uptake
            P_i_j=y_i_j*P_j;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(component_info,P_i_j,Nc).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(component_info,P_i_j).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=component_info(i).Fs(P_i_j').*rho;
                end
            end 
            % ode 
            sumdqdt=0;
            for i=1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));%solid phase mass balance
                sumdqdt=sumdqdt+dydt(Nf*(j-1)+Nc+i);
            end

             dydt(Nf*(j-1)+2*Nc+1)=(-1/delta_z)*(P_jph*v_jph-P_jmh*v_jmh)-((1-eps)/eps)*R*T*sumdqdt;% total mass balance

            for i = 1: Nc
                dydt(Nf*(j-1)+i)=(1/P_j)*Dax*(1/delta_z)*(P_jph*(y_i_jph(i)-y_i_j(i))/(delta_z/2)-P_jmh*(y_i_j(i)-y_i_jm1(i))/delta_z)...
                -(1/P_j)*(1/delta_z)*(y_i_jph(i)*P_jph*v_jph-y_i_jmh(i)*P_jmh*v_jmh)...
                -(R*T/P_j)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i)-(y_i_j(i)/P_j)*dydt(Nf*(j-1)+2*Nc+1);% gas phase mass balance
            end 
                
           
        else

            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            P_j=y(Nf*(j-1)+2*Nc+1);
            
            %j-2
            y_i_jm2=y(Nf*(j-3)+1:Nf*(j-3)+Nc);
            P_jm2=y(Nf*(j-3)+2*Nc+1);

            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);
            P_jm1=y(Nf*(j-2)+2*Nc+1);

            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);
            P_jp1=y(Nf*j+2*Nc+1);


            %j+0.5
            y_i_a0_j=(2./3)./((y_i_jp1-y_i_j+Delta).^4);
            y_i_a1_j=(1./3)/((y_i_j-y_i_jm1+Delta).^4);
            y_i_jph=(y_i_a0_j./(y_i_a0_j+y_i_a1_j)).*(y_i_j+y_i_jp1)./2+(y_i_a1_j./(y_i_a0_j+y_i_a1_j)).*(3*y_i_j./2-y_i_jm1./2);
            P_a0_j=(2./3)./((P_jp1-P_j+Delta).^4);
            P_a1_j=(1./3)./((P_j-P_jm1+Delta).^4);
            P_jph=(P_a0_j./(P_a0_j+P_a1_j)).*(P_j+P_jp1)./2+(P_a1_j./(P_a0_j+P_a1_j)).*(3*P_j./2-P_jm1./2);
            v_jph=(-4./150).*((eps./(1-eps)).^2)*(r_p.^2/mu).*(1./delta_z).*(P_jp1-P_j);

            %j-0.5
            y_i_a0_jm1=(2./3)./((y_i_j-y_i_jm1+Delta).^4);
            y_i_a1_jm1=(1./3)/((y_i_jm1-y_i_jm2+Delta).^4);
            y_i_jmh=(y_i_a0_jm1./(y_i_a0_jm1+y_i_a1_jm1)).*(y_i_jm1+y_i_j)./2+(y_i_a1_jm1./(y_i_a0_jm1+y_i_a1_jm1)).*(3*y_i_jm1./2-y_i_jm2./2);
            P_a0_jm1=(2./3)./((P_j-P_jm1+Delta).^4);
            P_a1_jm1=(1./3)./((P_jm1-P_jm2+Delta).^4);
            P_jmh=(P_a0_jm1./(P_a0_jm1+P_a1_jm1)).*(P_jm1+P_j)./2+(P_a1_jm1./(P_a0_jm1+P_a1_jm1)).*(3*P_jm1./2-P_jm2./2);
            v_jmh=(-4./150).*((eps./(1-eps)).^2).*(r_p.^2./mu)*(1./delta_z).*(P_j-P_jm1);

            % equilibrium uptake
            P_i_j=y_i_j*P_j;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(component_info,P_i_j,Nc).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(component_info,P_i_j).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=component_info(i).Fs(P_i_j').*rho;
                end
            end 
            % ode 
            sumdqdt=0;
            for i=1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));%solid phase mass balance
                sumdqdt=sumdqdt+dydt(Nf*(j-1)+Nc+i);
            end

             dydt(Nf*(j-1)+2*Nc+1)=(-1/delta_z)*(P_jph*v_jph-P_jmh*v_jmh)-((1-eps)/eps)*R*T*sumdqdt;% total mass balance

            for i = 1: Nc
                dydt(Nf*(j-1)+i)=(1/P_j)*Dax*(1/delta_z)*(P_jph*(y_i_jp1(i)-y_i_j(i))/delta_z-P_jmh*(y_i_j(i)-y_i_jm1(i))/delta_z)...
                -(1/P_j)*(1/delta_z)*(y_i_jph(i)*P_jph*v_jph-y_i_jmh(i)*P_jmh*v_jmh)...
                -(R*T/P_j)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i)-(y_i_j(i)/P_j)*dydt(Nf*(j-1)+2*Nc+1);% gas phase mass balance
            end 
                
           

        end 


    end 

end 
