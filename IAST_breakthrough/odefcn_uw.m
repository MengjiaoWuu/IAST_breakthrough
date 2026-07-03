function dydt=odefcn_uw(t,y,Nc,y_i,v_f,p_t,T,eps,k_i,Dax,rho,Nz,delta_z,isotherm_data_type,interpolation_method,varargin)

    Nf=2*Nc;% number of function
    c_f=MF_to_C(y_i,p_t,T);
    dydt=zeros(Nf*Nz,1);
    Delta=10^-(10);
    if strcmp(isotherm_data_type,'Isotherm Model')==1 && Nc==1 
            isotherm_type=varargin{1};
            b=varargin{2};
    elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1 && Nc==1 
            isotherm_data=varargin{1};
    else
            Fs=varargin{1};
    end


    for j = 1:Nz
        if j ==1  
            % j
            c_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
        
            %j+1
            c_i_jp1=y(Nf*j+1:Nf*j+Nc);

            %j+0.5
            c_i_jph=c_i_j;

            %j-0.5=1-0.5=0.5
%           c_i_jmh=(Dax.*c_i_j+v_f.*(delta_z./2).*c_f)./(v_f.*(delta_z./2)+Dax);
            c_i_jmh=c_f;
       
            % equilibrium uptake
            y_i_j=C_to_MF(c_i_j,p_t,T);
            P_i_j=y_i_j*p_t;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(isotherm_type{1},b{1},P_i_j,1).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(isotherm_data,P_i_j,interpolation_method).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=Fs{i}(P_i_j').*rho;
                end
            end 
            % ode 
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((c_i_jp1(i)-c_i_j(i))/delta_z-(c_i_j(i)-c_i_jmh(i))/(delta_z/2))...
                    -(v_f/delta_z)*(c_i_jph(i)-c_i_jmh(i))-((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        
        elseif j==Nz
            % j
            c_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);

        
            %j-1
            c_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);

            %j+1
            c_i_jp1=c_i_j;

            %j+0.5
            c_i_jph=c_i_j;

            %j-0.5
            c_i_jmh=c_i_jm1;
            
       
            % equilibrium uptake
            y_i_j=C_to_MF(c_i_j,p_t,T);
            P_i_j=y_i_j*p_t;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(isotherm_type{1},b{1},P_i_j,1).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(isotherm_data,P_i_j,interpolation_method).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=Fs{i}(P_i_j').*rho;
                end
            end 
            % ode 
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((c_i_jph(i)-c_i_j(i))/(delta_z/2)-(c_i_j(i)-c_i_jm1(i))/delta_z)...
                    -(v_f/delta_z)*(c_i_jph(i)-c_i_jmh(i))-((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        else
            % j
            c_i_j=y(Nf*(j-1)+1:Nf*(j-1)+Nc);
            q_i_j=y(Nf*(j-1)+Nc+1:Nf*(j-1)+2*Nc);
            

            %j-1
            c_i_jm1=y(Nf*(j-2)+1:Nf*(j-2)+Nc);

            %j+1
            c_i_jp1=y(Nf*j+1:Nf*j+Nc);

            
            %j+0.5
            c_i_jph=c_i_j;

            %j-0.5
            c_i_jmh=c_i_jm1;
            
       
            % equilibrium uptake
            y_i_j=C_to_MF(c_i_j,p_t,T);
            P_i_j=y_i_j*p_t;
            if Nc==1
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i_star_j=isotherm(isotherm_type{1},b{1},P_i_j,1).*rho;
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i_star_j=isotherm_interpolation(isotherm_data,P_i_j,interpolation_method).*rho;
                end 
            else
                for i = 1 : Nc
                    q_i_star_j(i)=Fs{i}(P_i_j').*rho;
                end
            end 
            % ode 
            for i = 1:Nc
                dydt(Nf*(j-1)+Nc+i)=k_i(i)*(q_i_star_j(i)-q_i_j(i));% solid phase mass balance
                dydt(Nf*(j-1)+i)=(Dax/delta_z)*((c_i_jp1(i)-c_i_j(i))/delta_z-(c_i_j(i)-c_i_jm1(i))/delta_z)...
                    -(v_f/delta_z)*(c_i_jph(i)-c_i_jmh(i))-((1-eps)/eps)*dydt(Nf*(j-1)+Nc+i);% gas phase mass balance
            end 
        end
    end 

end 