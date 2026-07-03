function[q_i,p_i_star]=IAST(isotherm_data_type,component_info,p_i,varargin)
    %This function uses the modified FastIAS algorithm to calculate gas mixture
    % uptake from gas composition and pure component isotherm model and
    % data 
    Nc=size(p_i,1);
    % isotherm points
    if strcmp(isotherm_data_type,'Isotherm Data Points')==1
        % error check
        for i = 1 :Nc
            if p_i > component_info(i).isotherm_data(end,1)
                    error(['Partial pressures of components in the gas ' ...
                        'mixture exceed the maximum pressure of their ' ...
                        'pure component isotherm data.']);
            elseif component_info(i).isotherm_data(1,1:2) ~=[0,0]
             component_info(i).isotherm_data=[[0,0];component_info(i).isotherm_data];
            end 
        end 
        % NQ_IAST
        [psi_bins_i,slope_i,intercept_i]=NQ_psi_bins(component_info,Nc);
        [q_i,p_i_star,q_i_star]=NQ_FastIAS(component_info,p_i,Nc,psi_bins_i,slope_i,intercept_i);
        % expolation check
        for i = 1:Nc 
            if component_info(i).isotherm_data(end,1)<p_i_star(i)
                warning(["During the numerical quadrature IAST calculation, "  +...
                    "the partial pressures of components in the theoretical"  +...
                    " adsorbed liquid phase exceed the maximum pressure of their" + ...
                    " pure component isotherm data. Constant slope extrapolation" + ...
                    " is used to calculate the uptake and reduced spreading pressure." + ...
                    " Please remeasure the pure component isotherm up to p_i_star" + ...
                    " to avoid extrapolation and ensure better accuracy."])
                p_i_star;
                break
            end 
    
        end 
        % figure check
        if isempty(varargin)==0 && strcmp(varargin,'figure off')==1 
            return
        else 
            figure
            newcolors = [    0    0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];
            for i = 1:Nc 
                plot (component_info(i).isotherm_data(:,1),component_info(i).isotherm_data(:,2),'o','MarkerSize',8,'LineWidth',1,'DisplayName',component_info(i).adsorbate,'Color',newcolors(i,:))
                hold on 
                xline(p_i(i),"-",'LineWidth',2.0,'DisplayName',append('p_{',component_info(i).adsorbate,'}'),'Color',newcolors(i,:))
                
                hold on
                xline(p_i_star(i),"--",'LineWidth',2.0,'DisplayName',append('p*_{',component_info(i).adsorbate,'}'),'Color',newcolors(i,:))
                hold on
            end 
        
        
            hold off
            xlabel('Pressure (Pa)')
            ylabel('Uptake (mol/kg)')
        
            ax = gca; 
            ax.FontSize = 15;
            ax.FontWeight='bold';
            hold off
            legend boxoff
            legend('Location','northeastoutside')
            legend show

        end 
    % isotherm model
    elseif strcmp(isotherm_data_type,'Isotherm Model')==1

        [q_i]=FastIAS(component_info,p_i,Nc);
        p_i_star=[];

    
    end 
end 