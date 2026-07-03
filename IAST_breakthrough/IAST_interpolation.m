function [component_info]=IAST_interpolation(isotherm_data_type,component_info)
    % This function generates spline interpolation functions of gas mixture
    % grids calculated with IAST.

    % check error 
    y_if=vertcat(component_info.y_if);
    Nc=size(component_info,1);
    if strcmp(isotherm_data_type,'Isotherm Data Points')==1
        for i = 1: Nc
            if  component_info(i).p_int_max> component_info(i).isotherm_data(end,1)
                error(['Maximum component partial pressures of the mixture' ...
                    ' isotherm grids exceed the maximum pressure of their' ...
                    ' pure component isotherm data.']);
            elseif component_info(i).isotherm_data(1,1:2) ~=[0,0]
                 component_info(i).isotherm_data=[[0,0];component_info(i).isotherm_data{i}];
            end 
        end
    end  
    % generate interpolation function 
    p_int_max=vertcat(component_info.p_int_max);
    interpolation_size=vertcat(component_info.interpolation_size);
    interpolation_method='spline';



    for i = 1: Nc
        n_level{i}=ceil(p_int_max(i)/interpolation_size(i));
        p_spread{i}=linspace(0,n_level{i}*interpolation_size(i),n_level{i}+1);
    end 
    grid_size=[n_level{:}]+1;
    [ps{1:Nc}] = ndgrid(p_spread{1:Nc});
    
    out = cell2mat(cellfun(@(x)x(:),ps,'un',0));
    disp(['The IAST calculated mixture isotherm grid size is ',num2str(grid_size)])

    if strcmp(isotherm_data_type,'Isotherm Data Points')==1
        [psi_bins_i,slope_i,intercept_i]=NQ_psi_bins(component_info,Nc);
    end 
    p_i_star=zeros(Nc,1);
    for i =1 : size(out,1)
            p=out(i,:)';
            if sum(p)==0
                q_i=p;
            elseif nnz(p)==1
                p_index=find(p~=0);
                q_i=p;
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    q_i(p_index)=isotherm(component_info(p_index),p(p_index),1);
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    q_i(p_index)=isotherm_interpolation(component_info(p_index),p(p_index));
                end
            elseif nnz(p)==Nc
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    [q_i]=FastIAS(component_info,p,Nc);
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    [q_i,pp_i_star]=NQ_FastIAS(component_info,p,Nc,psi_bins_i,slope_i,intercept_i);
                    p_i_star=max([p_i_star,pp_i_star],[],2);
                end 
            elseif nnz(p)<Nc
                p_index=find(p~=0);
                q_i=zeros(1,Nc);
                pp_i_star=zeros(1,Nc);
                new_p=p(p_index);
                new_Nc=size(new_p,1);
                new_component_info=component_info(p_index);
                if strcmp(isotherm_data_type,'Isotherm Model')==1
                    [new_q_i]=FastIAS(new_component_info,new_p,new_Nc);
                elseif strcmp(isotherm_data_type,'Isotherm Data Points')==1
                    [new_psi_bins_i,new_slope_i,new_intercept_i]=NQ_psi_bins(new_component_info,new_Nc);
                    [new_q_i,new_pp_i_star]=NQ_FastIAS(component_info,p,new_Nc,new_psi_bins_i,new_slope_i,new_intercept_i);
                    pp_i_star(p_index)=new_pp_i_star;
                    p_i_star=max([p_i_star,pp_i_star],[],2);
                end 
                q_i(p_index)=new_q_i;

            end 
            uptake_M(i,:)=q_i;
            if rem(i,10000)==0
                disp(['Finished points:',num2str(i),' Total points:',num2str(size(out,1))])
            end 

    end 
    disp('IAST calculation is done')
    for i = 1: Nc
        F = griddedInterpolant(ps{1:Nc},reshape(uptake_M(:,i),grid_size),interpolation_method);
        Fs{i}=F;
        component_info(i).Fs=F;
    end 
    % check expolation 
    if strcmp(isotherm_data_type,'Isotherm Data Points')==1
        for i=1: Nc
        if component_info(i).isotherm_data(end,1)<p_i_star(i)
            warning(["During the numerical quadrature IAST calculation, "  +...
                "the partial pressures of components in the theoretical"  +...
                " adsorbed liquid phase exceed the maximum pressure of their" + ...
                " pure component isotherm data. Constant slope extrapolation" + ...
                " is used to calculate the uptake and reduced spreading pressure." + ...
                " Please remeasure the pure component isotherm up to p_i_star" + ...
                " to avoid extrapolation and ensure better accuracy."])
            p_i_star
            break
        end 

        end 
    end 
end 