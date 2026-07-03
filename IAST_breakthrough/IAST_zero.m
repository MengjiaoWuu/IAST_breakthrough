function [q_i,p_i_star]=IAST_zero(isotherm_data_type,component_info,p)
    % this function accept gas composition with some zero partial pressure 
    Nc=size(p,1);

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
        p_i_star=[];
    elseif nnz(p)==Nc
        [q_i,p_i_star]=IAST(isotherm_data_type,component_info,p,'figure off');
    elseif nnz(p)<Nc
        p_index=find(p~=0);
        q_i=zeros(1,Nc);
        p_i_star=zeros(1,Nc);
        new_p=p(p_index);
        new_component_info=component_info(p_index);
        [new_q_i,new_p_i_star]=IAST(isotherm_data_type, new_component_info,new_p,'figure off');
        if strcmp(isotherm_data_type,'Isotherm Model')==1
            p_i_star=new_p_i_star;
        else 
            p_i_star(p_index)=new_p_i_star;
        end 
        q_i(p_index)=new_q_i;
    end 

end 