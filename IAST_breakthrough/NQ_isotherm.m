function [q_i]=NQ_isotherm(component_info,slope_i,intercept_i,p_i,Nc)
    % calculate uptake using linear interpolation 
    for i =1: Nc
        p=p_i(i);
        p_data=component_info(i).isotherm_data(:,1);
        
        slope=slope_i{i};
        intercept=intercept_i{i};
        for j =2: size(p_data,1)
            if p <p_data(1)
                q=0;
                break
            elseif p <= p_data(j)
                q=intercept(j-1)+p.*slope(j-1);
                break
             elseif p > p_data(end)
                % % constant slope
                q=intercept(end)+p.*slope(end);

                % end
                break     
            end 
        end 
        q_i(i,1)=q;
        
    end 
end 

