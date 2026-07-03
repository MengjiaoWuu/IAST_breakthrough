function [c_i]=MF_to_C(y_i,p_t,T)
    %This function converts molar ratio to concentration 
    %P:total pressure Pa
    %T:temperature K
    R_g=8.314;% m^3 Pa/(mol K) gas constant 
    c_t=p_t./(R_g.*T);% total concentration mol/m^3
    c_i=y_i.*c_t;% mol/m^3
end 