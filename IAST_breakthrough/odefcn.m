function dydt=odefcn(t,y,v,P,T,eps,rho,Dax,Nz,delta_z,isotherm_data_type,component_info)

    % parameters
    Nc=size(component_info,1);
    y_if=vertcat(component_info.y_if);
    k_i=vertcat(component_info.k_i);
    R=8.3144598;% gas constant (m^3 Pa/(mol k))


    interpolation_method='spline';
    Nf=2*Nc;% number of function
    dydt=zeros(Nf*Nz,1);
    Delta=10^-(10);



    for j = 1:Nz
        if j ==1  
            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
        
            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);

            %j+0.5
            y_i_jph=y_i_j;

            %j-0.5=1-0.5=0.5
            y_i_jmh=(Dax.*y_i_j+v.*(delta_z./2).*y_if)./(v.*(delta_z./2)+Dax);

            % equilibrium uptake
            P_i_j=y_i_j*P;
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
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((y_i_jp1(i)-y_i_j(i))/delta_z-(y_i_j(i)-y_i_jmh(i))/(delta_z/2))...
                    -(v/delta_z)*(y_i_jph(i)-y_i_jmh(i))-(R*T/P)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        elseif j ==2
           % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            
            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);

            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);
            
            %j+0.5
            y_i_jph=y_i_j;

            %j-0.5=1-0.5=0.5
            y_i_jmh=y_i_jm1;
       
            % equilibrium uptake
            P_i_j=y_i_j*P;
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
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((y_i_jp1(i)-y_i_j(i))/delta_z-(y_i_j(i)-y_i_jm1(i))/delta_z)...
                    -(v/delta_z)*(y_i_jph(i)-y_i_jmh(i))-(R*T/P)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 


        elseif j==Nz
            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);

            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);

            %j+0.5=Nz+0.5(end)
            y_i_jph=y_i_j;
            %j-0.5
            y_i_jmh=y_i_jm1;
       
            % equilibrium uptake
            P_i_j=y_i_j*P;
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
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((y_i_jph(i)-y_i_j(i))/(delta_z/2)-(y_i_j(i)-y_i_jm1(i))/delta_z)...
                    -(v/delta_z)*(y_i_jph(i)-y_i_jmh(i))-(R*T/P)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        else
            % j
            y_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            
            %j-2
            y_i_jm2=y(Nf*(j-3)+1:Nf*(j-3)+Nc);

            %j-1
            y_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);

            %j+1
            y_i_jp1=y(Nf*j+1:Nf*j+Nc);

            
            %j+0.5
            y_i_a0_j=(2/3)/((y_i_jp1-y_i_j+Delta).^4);
            y_i_a1_j=(1/3)/((y_i_j-y_i_jm1+Delta).^4);
            y_i_jph=(y_i_a0_j/(y_i_a0_j+y_i_a1_j))*(y_i_j+y_i_jp1)/2+(y_i_a1_j/(y_i_a0_j+y_i_a1_j))*(3*y_i_j/2-y_i_jm1/2);

            %j-0.5
            y_i_a0_jm1=(2/3)/((y_i_j-y_i_jm1+Delta).^4);
            y_i_a1_jm1=(1/3)/((y_i_jm1-y_i_jm2+Delta).^4);
            y_i_jmh=(y_i_a0_jm1/(y_i_a0_jm1+y_i_a1_jm1))*(y_i_jm1+y_i_j)/2+(y_i_a1_jm1/(y_i_a0_jm1+y_i_a1_jm1))*(3*y_i_jm1/2-y_i_jm2/2);
            
       
            % equilibrium uptake
            P_i_j=y_i_j*P;
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
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((y_i_jp1(i)-y_i_j(i))/delta_z-(y_i_j(i)-y_i_jm1(i))/delta_z)...
                    -(v/delta_z)*(y_i_jph(i)-y_i_jmh(i))-(R*T/P)*((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        end
    end 

end 