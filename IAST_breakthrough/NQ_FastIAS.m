function[q_i,p_i_star,q_i_star]=NQ_FastIAS(component_info,p_i,Nc,psi_bins_i,slope_i,intercept_i)
    % This function calculates mixture uptake with pure component isotherm
    % data 
    % inital guess
    p_t=sum(p_i);
    y_i=p_i./p_t;


    for i=1:Nc
        H_i(i,1)=slope_i{i}(1);
    end 
    H_ave=sum(y_i.*H_i);
    col(1,:)=p_t.*H_ave./H_i;
    col(2,:)=p_t;
    p_i_star=min(col)';
  
    psi_i=NQ_RGP(component_info,psi_bins_i,slope_i,intercept_i,p_i_star,Nc);
    
    % iteration 
    test=1;
    while test >= 10^(-13)
        % intial matrix
        phi_i= NQ_isotherm(component_info,slope_i,intercept_i,p_i_star,Nc)./p_i_star;
        M_phi=diag(phi_i);
        M_phi(:,end)=-phi_i(end);
        M_phi(end,:)=-y_i.*p_t./(p_i_star.^2);       
        G=zeros(Nc,1);
        for i = 1 : Nc
            
            if i~=Nc
                G(i)=psi_i(i)-psi_i(Nc);
            else
                G(i)=sum((y_i.*p_t)./p_i_star)-1;
            end
        end
        
        % update matrix
        M_phi(Nc,Nc)=M_phi(Nc,Nc)-M_phi(Nc,1:Nc-1)./diag(M_phi(1:Nc-1,1:Nc-1))'*(M_phi(1:Nc-1,Nc));
        G(Nc)=G(Nc)-M_phi(Nc,1:Nc-1)./diag(M_phi(1:Nc-1,1:Nc-1))'*(G(1:Nc-1));
        
        % delta
        delta=zeros(Nc,1);
        
        for i= Nc:-1:1
            if i==Nc
                delta(i)=G(i)/M_phi(i,i);
            else
                delta(i)=(G(i)-delta(Nc)*M_phi(i,Nc))/(M_phi(i,i));
            end
        end
        
        new_p_i_star=p_i_star-delta;
    
        for i =1: Nc
            if new_p_i_star(i) < 0
                new_p_i_star(i)=p_i_star(i)/2;
            end
        end
        % max_p_i_star=max([p_i_star,new_p_i_star],[],2);
        p_i_star=new_p_i_star;
        psi_i=NQ_RGP(component_info,psi_bins_i,slope_i,intercept_i,p_i_star,Nc);
        test=sum(abs(delta./p_i_star));
        
    end
    
    x_i_final=y_i.*p_t./p_i_star;
    psi_i_final=NQ_RGP(component_info,psi_bins_i,slope_i,intercept_i,p_i_star,Nc);
    q_i_star=NQ_isotherm(component_info,slope_i,intercept_i,p_i_star,Nc);
    q_t=1/sum(x_i_final./q_i_star);
    q_i=x_i_final.*q_t;
end

