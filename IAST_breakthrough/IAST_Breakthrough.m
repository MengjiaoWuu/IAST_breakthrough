function [output]=IAST_Breakthrough(adsorbent,isotherm_data_type,component_info,time,Nz,v,P,L,T,eps,rho,Dax,model_type,show_figure,varargin)
    % This function simulates breakthrough curves with IAST using  pregenerated interpolation function 
    adsorbate=vertcat(component_info.adsorbate);
    y_if=vertcat(component_info.y_if);
    Nc=size(component_info,1);
    delta_z=L/Nz;
    zbd=linspace(0,L,Nz+1);
    zs=zbd(1:Nz)+delta_z/2;
    c_if=MF_to_C(y_if,P,T);
    p_if=y_if.*P;
    tspan=[0,time];

    if strcmp(model_type,'Pv_constant')==1
        y0=zeros(1,2*Nc*Nz); % initial bed condition (start with a  clean bed)
        %ode
         [t,y] =ode15s(@(t,y) odefcn(t,y,v,P,T,eps,rho,Dax,Nz,delta_z,isotherm_data_type,component_info), tspan,y0);
        % process result 
        for i=1:Nc
            y_i{i}=y(:,i:(2*Nc):(2*Nc)*(Nz-1)+i);
            q_i{i}=y(:,Nc+i:(2*Nc):(2*Nc)*(Nz-1)+Nc+i)/rho;
            c_i{i}=MF_to_C(y_i{i},P,T);
            p_i{i}=y_i{i}.*P;
            breakthrough_curve_c{i}(:,1)=t;
            breakthrough_curve_c{i}(:,2)=c_i{i}(:,end)/c_if(i);
            breakthrough_curve_p{i}(:,1)=t;
            breakthrough_curve_p{i}(:,2)=p_i{i}(:,end)/p_if(i);
        end 
        

        if strcmp(show_figure,'figure on')==1 

            figure 
            for i=1:Nc
                plot (t./60,p_i{i}(:,end)/p_if(i),'LineWidth',2.0)
                hold on
            end
            xlabel('Time(min)')
            ylabel('p_i/p_{ifeed}')
            ax = gca; 
            ax.FontSize = 15;
            xlim([0 tspan(2)./60])
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off 
    
    
            figure 
            for i=1:Nc
                plot (t./60,c_i{i}(:,end)/c_if(i),'LineWidth',2.0)
                hold on
            end
            xlabel('Time(min)')
            ylabel('c_i/c_{ifeed}')
            xlim([0 tspan(2)./60])
            ax = gca; 
            ax.FontSize = 15;
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off 
    
    
            
            figure 
            for i=1:Nc
                plot (zs,reshape(q_i{i}(end,:),[1,size(zs,2)]),'LineWidth',2.0)
                hold on
            end
            xlabel('Bed Length(m)')
            ylabel('Amount Adsorbed(mol/kg)')
    
            ax = gca; 
            ax.FontSize = 15;
            xlim([0 L])
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off 

        end 
        unit=struct('adsorbate',[]','breakthrough_curve_pioverpf',[],'breakthrough_curve_ciovercf',[],'p_i',"Pa",'c_i',"mol/m^3",'q_i',"mol/kg",'t',"s",'zs',"m");
        result=struct('adsorbate',cellstr(adsorbate)','breakthrough_curve_pioverpf',breakthrough_curve_p,'breakthrough_curve_ciovercf',breakthrough_curve_c,'p_i',p_i,'c_i',c_i,'q_i',q_i,'t',t,'zs',zs);
        output=[unit,result];

    elseif strcmp(model_type,'Pv_drop')==1
        r_p=varargin{1};
        mu=varargin{2};
        y0=zeros(1,(1+2*Nc)*Nz); % initial bed condtion (start with a  clean bed) 
        y0(2*Nc+1:2*Nc+1:(1+2*Nc)*Nz)=P;% initial bed condtion (start with a  clean bed with regulated pressure)


        %ode
        [t,y] =ode15s(@(t,y) odefcn_PD(t,y,v,P,T,eps,rho,Dax,Nz,delta_z,isotherm_data_type,component_info,r_p,mu), tspan,y0);
        p=y(:,(2*Nc+1):(2*Nc+1):(2*Nc+1)*Nz);
        for i=1:Nc
            y_i{i}=y(:,i:(2*Nc+1):(2*Nc+1)*(Nz-1)+i);
            q_i{i}=y(:,Nc+i:(2*Nc+1):(2*Nc+1)*(Nz-1)+Nc+i)/rho;
            c_i{i}=MF_to_C(y_i{i},p,T);
            p_i{i}=y_i{i}.*p;
            breakthrough_curve_c{i}(:,1)=t;
            breakthrough_curve_c{i}(:,2)=c_i{i}(:,end)/c_if(i);
            breakthrough_curve_p{i}(:,1)=t;
            breakthrough_curve_p{i}(:,2)=p_i{i}(:,end)/p_if(i);

        end   
     

        if strcmp(show_figure,'figure on')==1
            figure 
            for i=1:Nc
                plot (t./60,p_i{i}(:,end)/p_if(i),'LineWidth',2.0)
                hold on
            end
            xlabel('Time(min)')
            ylabel('p_i/p_{ifeed}')
            ax = gca; 
            ax.FontSize = 15;
            xlim([0 tspan(2)./60])
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off 
    
            figure 
            for i=1:Nc
                plot (t./60,c_i{i}(:,end)/c_if(i),'LineWidth',2.0)
                hold on
            end
            xlabel('Time(min)')
            ylabel('c_i/c_{ifeed}')
            ax = gca; 
            ax.FontSize = 15;
            xlim([0 tspan(2)./60])
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off 
            
    
            figure
            plot(zs,p(end,:),'LineWidth',2.0)
            ax = gca; 
            ax.FontSize = 15;
            xlabel('Bed Length(m)')
            ylabel('P(Pa)')
            xlim([0 L])
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off
            
            figure 
            for i=1:Nc
                plot (zs',q_i{i}(end,:),'LineWidth',2.0)
                hold on
            end
            ax = gca; 
            ax.FontSize = 15;
            xlabel('Bed Length(m)')
            ylabel('Amount Adsorbed(mol/kg)')
            xlim([0 L])
            legend(adsorbate)
            title(adsorbent,['Temperature:',num2str(T),' K  ','Pressure:',num2str(P),' Pa'])
            hold off
        end 

        unit=struct('adsorbate',[]','breakthrough_curve_pioverpf',[],'breakthrough_curve_ciovercf',[],'p_i',"Pa",'c_i',"mol/m^3",'q_i',"mol/kg",'t',"s",'zs',"m",'P',"Pa",'y_i',[]);
        result=struct('adsorbate',cellstr(adsorbate)','breakthrough_curve_pioverpf',breakthrough_curve_p,'breakthrough_curve_ciovercf',breakthrough_curve_c,'p_i',p_i,'c_i',c_i,'q_i',q_i,'t',t,'zs',zs,'P',p,'y_i',y_i);
        output=[unit,result];
    end 

end 
   












