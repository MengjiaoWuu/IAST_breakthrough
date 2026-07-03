function [psi_i]=NQ_RGP(component_info,psi_bins_i,slope_i,intercept_i,p_i,Nc)
    % calculate the reduced grand potential for a certain pressure with
    % pure isotherm data 
    for i =1: Nc
        p=p_i(i);
        p_data=component_info(i).isotherm_data(:,1);
        slope=slope_i{i};
        intercept=intercept_i{i};
        psi_bins=psi_bins_i{i};
        for j =1: size(p_data,1)
            if abs(p-p_data(j))<10^(-10)
                psi=sum(psi_bins(1:j-1));
                break
            elseif p < p_data(2)
                psi=psi_bins(1)*p/p_data(2);
                break
            elseif p <= p_data(j)
                psi=sum(psi_bins(1:j-2))+slope(j-1).*(p-p_data(j-1))+intercept(j-1).*log(p./p_data(j-1));
                break
            elseif p > p_data(end)
                last_point=size(p_data,1);
                %constant slope
                psi=sum(psi_bins(1:last_point-1))+slope(last_point-1).*(p-p_data(last_point))+intercept(last_point-1).*log(p./p_data(last_point));

                break
            end 
        end 
        psi_i(i,1)=psi;
        
    end 
end 

