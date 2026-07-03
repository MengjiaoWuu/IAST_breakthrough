function [psi]=RGP(component_info,p,Nc)
    % calculate the reduced grand potential for a certain pressure with
    % pure isotherm model 
    psi=zeros(Nc,1);
    for i =1 : Nc
        B=component_info(i).parameters;
        iso=component_info(i).isotherm_type;
        if isfield(component_info(i),'number_of_site')==1
            site=component_info(i).number_of_site;
        end 
        if strcmp(iso,'Langmuir')==1
            psi_i=0;   
            for j = 1:site 
                psi_i=psi_i+B(2*(j-1)+1)*log(1+B(2*(j-1)+2).*p(i));
            end 
        elseif strcmp(iso,'BET')==1
            psi_i=B(1)*log(1+(1-B(3)*p(i)+B(2)*p(i)/(1-B(3)*p(i))));
        elseif strcmp(iso,'Anti-Langmuir')==1
            psi_i=(-B(1)/B(2))*log(1-B(2)*p(i));
        elseif strcmp(iso,'Henry')==1
            psi_i=B(1)*p(i);
        elseif strcmp(iso,'Freundlich')==1
            psi_i=B(1)*B(2)*p(i)^(1/B(2));
        elseif strcmp(iso,'Sips')==1
            psi_i=0;
            for j = 1:site 
                psi_i=psi_i+B(3*(j-1)+1)*B(3*(j-1)+3)*log(1+(B(3*(j-1)+2)*p(i))^(1/B(3*(j-1)+3)));
            end 
        elseif strcmp(iso,'Langmuir-Freundlich')==1
            
            psi_i=0;
            for j = 1:site 
                psi_i=psi_i+(B(3*(j-1)+1)/B(3*(j-1)+3))*log(1+B(3*(j-1)+2)*p(i)^B(3*(j-1)+3));
            end 
        elseif strcmp(iso,'OBrien&Myers')==1
            psi_i=B(1)*(log(1+B(2)*p(i))+(B(3)^2)*B(2)*p(i)/(2*(1+B(2)*p(i))^2));
        elseif strcmp(iso,'Quadratic')==1
            psi_i=B(1)*log(1+B(2)*p(i)+B(3)*p(i)^2);
        elseif strcmp(iso,'Asymptotic Temkin')==1
            psi_i=B(1)*(log(1+B(2)*p(i))-(1/2)*B(3)*(B(2)*p(i)/(1+B(2)*p(i)))^2);
        end
        psi(i,1)=psi_i;
    end
end
