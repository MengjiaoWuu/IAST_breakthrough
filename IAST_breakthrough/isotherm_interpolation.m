function q=isotherm_interpolation(component_info,p)
% This function calculates single gas uptake form pure isotherm data points with spline interpolation 
    isotherm_data=component_info(1).isotherm_data;

    x=isotherm_data(:,1);
    v=isotherm_data(:,2);
    q = interp1(x,v,p,'spline');
end 