function [psi_bins_i,slope_i,intercept_i]=NQ_psi_bins(component_info,Nc)
    % find the slope and intercept of linear interpolation 
    % find the reduced grand potential between each two data points

    if Nc==1
        p=component_info(1).isotherm_data(:,1);
        q=component_info(1).isotherm_data(:,2);
        p_diff=diff(p);
        q_diff=diff(q);
        slope=q_diff./p_diff;
        intercept=q(2:end)-p(2:end).*slope;
        slope_i{1}=slope;
        intercept_i{1}=intercept;
        psi_bins_i=[];
    else 
        for i =1 :Nc
            p=component_info(i).isotherm_data(:,1);
            q=component_info(i).isotherm_data(:,2);
            p_diff=diff(p);
            q_diff=diff(q);
            slope=q_diff./p_diff;
            intercept=q(2:end)-p(2:end).*slope;

            psi_bins=slope.*p_diff+intercept.*log(p(2:end)./p(1:end-1));
            psi_bins(1)=q(2);
            psi_bins_i{i}=psi_bins;
 
            slope_i{i}=slope;
            intercept_i{i}=intercept;
        end 
    end 
end



