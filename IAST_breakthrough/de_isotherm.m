function [phi]= de_isotherm(component_info,p,Nc)
    q=isotherm(component_info,p,Nc);
    phi=q./p;
end 