function[y_i]=C_to_MF(c_i,p_t,T)
    %This function converts concentration to molar fraction.
    %c_i:concentration mol/m^3
    %p_t:total pressure pa
    %T:temperature k
    R_g=8.314;% m^3 Pa/(mol k) gas constant 
    c_t=p_t./(R_g.*T);
    y_i=c_i./c_t;

end 