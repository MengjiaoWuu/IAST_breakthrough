function [q]= isotherm(component_info,p,Nc)
    % This function calculates uptake from pure component isotherm model
    if Nc==1
        p_length=size(p,2);
        B=component_info.parameters;
        iso=component_info.isotherm_type;
        if isfield(component_info,'number_of_site')==1
            site=component_info.number_of_site;
        end 
        q=zeros(1,p_length);
        if strcmp(iso,'Langmuir')==1
            for j = 1:site 
                q=q+B(2.*(j-1)+1).*B(2.*(j-1)+2).*p./(1+B(2*(j-1)+2).*p);
            end 
        elseif strcmp(iso,'BET')==1
            q=B(1).*B(2).*p./((1-B(3).*p).*(1-B(3).*p+B(2).*p));
        elseif strcmp(iso,'Anti-Langmuir')==1
            q=B(1)*p./(1-B(2).*p);
        elseif strcmp(iso,'Henry')==1
            q=B(1).*p;
        elseif strcmp(iso,'Freundlich')==1
            q=B(1).*p.^(1./B(2));
        elseif strcmp(iso,'Sips')==1
            q=zeros(1,p_length);
            for j = 1:site 
                q=q+B(3.*(j-1)+1).*(B(3.*(j-1)+2).*p).^(1/B(3.*(j-1)+3))./(1+(B(3.*(j-1)+2).*p).^(1./B(3.*(j-1)+3)));
            end 
        elseif strcmp(iso,'Langmuir-Freundlich')==1
            q=zeros(1,p_length);
            for j = 1:site 
                q=q+B(3.*(j-1)+1).*B(3.*(j-1)+2).*p.^B(3.*(j-1)+3)./(1+B(3.*(j-1)+2).*p.^B(3.*(j-1)+3));
            end 
        elseif strcmp(iso,'OBrien&Myers')==1
            q=B(1).*(B(2).*p./(1+B(2).*p)+(B(3).^2).*B(2).*p.*(1-B(2).*p)/(2.*(1+B(2).*p).^3));
        elseif strcmp(iso,'Quadratic')==1
            q=B(1).*(B(2).*p+2.*B(3).*p.^2)./(1+B(2).*p+B(3).*p.^2);
        elseif strcmp(iso,'Asymptotic Temkin')==1
            q=B(1).*B(2).*p./(1+B(2).*p)+(B(1).*B(3).*(B(2).*p./(1+B(2).*p)).^2).*(B(2).*p./(1+B(2).*p)-1);
        end
    else 
        q=zeros(Nc,1);
        for i =1 : Nc
            B=component_info(i).parameters;
            iso=component_info(i).isotherm_type;
            if isfield(component_info(i),'number_of_site')==1
                site=component_info(i).number_of_site;
            end 
            if strcmp(iso,'Langmuir')==1
                q_i=0;   
                for j = 1:site 
                    q_i=q_i+B(2.*(j-1)+1).*B(2.*(j-1)+2).*p(i)./(1+B(2*(j-1)+2).*p(i));
                end 
            elseif strcmp(iso,'BET')==1
                q_i=B(1)*B(2)*p(i)/((1-B(3)*p(i))*(1-B(3)*p(i)+B(2)*p(i)));
            elseif strcmp(iso,'Anti-Langmuir')==1
                q_i=B(1)*p(i)/(1-B(2)*p(i));
            elseif strcmp(iso,'Henry')==1
                q_i=B(1)*p(i);
            elseif strcmp(iso,'Freundlich')==1
                q_i=B(1).*p(i).^(1/B(2));
            elseif strcmp(iso,'Sips')==1
                q_i=0;
                for j = 1:site 
                    q_i=q_i+B(3*(j-1)+1)*(B(3*(j-1)+2)*p(i))^(1/B(3*(j-1)+3))/(1+(B(3*(j-1)+2)*p(i))^(1/B(3*(j-1)+3)));
                end 
            elseif strcmp(iso,'Langmuir-Freundlich')==1
                q_i=0;
                for j = 1:site 
                    q_i=q_i+B(3*(j-1)+1)*B(3*(j-1)+2)*p(i)^B(3*(j-1)+3)/(1+B(3*(j-1)+2)*p(i)^B(3*(j-1)+3));
                end 
            elseif strcmp(iso,'OBrien&Myers')==1
                q_i=B(1)*(B(2)*p(i)/(1+B(2)*p(i))+(B(3)^2)*B(2)*p(i)*(1-B(2)*p(i))/(2*(1+B(2)*p(i))^3));
            elseif strcmp(iso,'Quadratic')==1
                q_i=B(1).*(B(2).*p(i)+2.*B(3).*p(i).^2)/(1+B(2).*p(i)+B(3).*p(i).^2);
            elseif strcmp(iso,'Asymptotic Temkin')==1
                q_i=B(1)*B(2)*p(i)/(1+B(2)*p(i))+B(1)*B(3)*(B(2)*p(i)/(1+B(2)*p(i)))^2*(B(2)*p(i)/(1+B(2)*p(i))-1);
            end
            q(i,1)=q_i;
        end
    end
end
