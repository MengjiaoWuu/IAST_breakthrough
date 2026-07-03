function [E]=IAST_interpolation_accuracy(P,isotherm_data_type,component_info)
    Np=500;
    
    y_if=vertcat(component_info.y_if);
    Nc=size(y_if,1);
    for i =1:Nc
        p_test(i,1:Np)=rand(1,Np).*y_if(i).*P;
    end
    
    if strcmp(isotherm_data_type,'Isotherm Data Points')==1
        for i =1:Nc
            isotherm_data{i}=component_info(i).isotherm_data;
        end
    end 
    for i =1:Nc
        for j=1:Np
            predict(i,j)=component_info(i).Fs(p_test(:,j)');
            if strcmp(isotherm_data_type,'Isotherm Data Points')==1
                [psi_bins_i,slope_i,intercept_i]=NQ_psi_bins(isotherm_data,Nc);
                [actual(:,j),p_i_star,q_i_star]=NQ_FastIAS(isotherm_data,p_test(:,j),Nc,psi_bins_i,slope_i,intercept_i);
            elseif strcmp(isotherm_data_type,'Isotherm Model')==1
                [actual(:,j)]=FastIAS(component_info,p_test(:,j),Nc);
            end 
        end 
    end 
    
    for i = 1:Nc 
          E(i,1)=100*mean(abs((actual(i,:)-predict(i,:)))./actual(i,:));
    end



end 