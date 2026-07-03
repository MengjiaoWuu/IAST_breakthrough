function[H]=Henry(Nc,component_info)
    % This function calculates Henry's constant for both isotherm model 


    p=0.01;% pa
    for i =1 : Nc
        B=component_info(i).parameters;
        iso=component_info(i).isotherm_type;
        site=component_info(i).number_of_site;
        if strcmp(iso(1),'Langmuir')==1  
            H_i=0;   
            for j = 1:site 
                H_i=H_i+B(2.*(j-1)+1).*B(2.*(j-1)+2);
            end 
        elseif strcmp(iso(1),'BET')==1
            H_i=B(1)*B(2)/(1-B(3));
        elseif strcmp(iso(1),'Anti-Langmuir')==1
            H_i=B(1);
        elseif strcmp(iso(1),'Henry')==1
            H_i=B(1);
        elseif strcmp(iso(1),'Freundlich')==1
            H_i=B(1).*p.^(1/B(2))/p;
        elseif strcmp(iso(1),'Sips')==1
            H_i=0;
            for j = 1:site 
                H_i=H_i+B(3*(j-1)+1)*(B(3*(j-1)+2)*p)^(1/B(3*(j-1)+3))/(1+(B(3*(j-1)+2)*p)^(1/B(3*(j-1)+3)))/p;
            end 
        elseif strcmp(iso(1),'Langmuir-Freundlich')==1
            
            H_i=0;
            for j = 1:site 
                H_i=H_i+B(3*(j-1)+1)*B(3*(j-1)+2)*p^B(3*(j-1)+3)/(1+B(3*(j-1)+2)*p^B(3*(j-1)+3))/p;
            end 
        elseif strcmp(iso(1),'OBrien&Myers')==1
            H_i=B(1)*B(2)*(1+(B(2)^2)/2);
        elseif strcmp(iso(1),'Quadratic')==1
            H_i=B(1)*B(2);
        elseif strcmp(iso(1),'Asymptotic Temkin')==1
            H_i=B(1)*B(2);
        end
        H(i,1)=H_i;
    end



end 