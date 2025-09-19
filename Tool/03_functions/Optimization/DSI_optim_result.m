function info_new = DSI_optim_result(var,info)

    switch info.version
        case 1
            net_var=4;
        case 2
            net_var=4;
    end

    
    switch info.type
        case 'Grid'
            var_grid = var([1:net_var]);
            bool_reload_vsc=false;
        case 'Electric'
            var_grid = var([1:net_var]);
            var_vsc = zeros(1,18);
            var_vsc([1:6]) =var([net_var+1:net_var+6]);
            bool_reload_vsc=true;
        case 'All'
            var_grid = var([1:net_var]);
            var_vsc = zeros(1,18);
            var_vsc([1:15]) =var([net_var+1:15+net_var]);
            bool_reload_vsc=true;
    end

    T_NET           = readtable(info.newmainfile,'Sheet','AC-NET');
    T_load          = readtable(info.newmainfile,'Sheet','load');   
    T_VSC           = readtable(info.newmainfile,'Sheet','VSC');  

    switch info.version
        case 1

            alpha_ii1=var_grid(1);
            alpha_ii2=var_grid(2);
        
            S_load=alpha_ii1*real(info.network.S_iny_base)+alpha_ii2*imag(info.network.S_iny_base)*1i;
            S_losses=info.network.S_iny_base-S_load;
        
            beta_jj1=var_grid(3);
            beta_jj2=var_grid(4);
        
            S_L1=beta_jj1*real(S_load)+beta_jj2*imag(S_load)*1i;
            S_L2=(1-beta_jj1)*real(S_load)+(1-beta_jj2)*imag(S_load)*1i;
            
            S_12=info.network.S_net_base-S_L1;
            
            I_12=conj(S_12/info.network.V_2_base);
            Z_eq=S_losses/abs(I_12)^2;
        
            T_load.P(1)=real(S_L1);
            T_load.Q(1)=imag(S_L1);
            T_load.P(2)=real(S_L2);
            T_load.Q(2)=imag(S_L2);
        
            T_NET.R(1)=real(Z_eq);
            T_NET.X(1)=imag(Z_eq);

        case 2

            % alpha_p=var_grid(1);
            % alpha_q=var_grid(2);
            % 
            % beta_1p=0;            
            % beta_2p=var_grid(3);
            % 
            % beta_1q=0;
            % beta_2q=var_grid(4);
            % 
            % gamma_p=var_grid(5);
            % gamma_q=var_grid(6);
            
            % 
            % S_load=alpha_p*real(info.network.S_iny_base)+alpha_q*imag(info.network.S_iny_base)*1i;
            % S_losses=info.network.S_iny_base-S_load;
            % 
            % S_L1=beta_1p*real(S_load)+beta_1q*imag(S_load)*1i;
            % S_L2=beta_2p*real(S_load)+beta_2q*imag(S_load)*1i;
            % S_L3=(1-beta_1p-beta_2p)*real(S_load)+(1-beta_1q-beta_2q)*imag(S_load)*1i;
            % 
            % S_losses_12=gamma_p*real(S_losses)+gamma_q*imag(S_losses)*1i;
            % S_losses_23=(1-gamma_p)*real(S_losses)+(1-gamma_q)*imag(S_losses)*1i;
            % 
            % S_12=info.network.S_net_base-S_L1;          
            % I_12=conj(S_12/info.network.V_2_base);
            % 
            % Z_eq_12=S_losses_12/abs(I_12)^2;
            % 
            % V_3_base=-I_12*Z_eq_12+info.network.V_2_base;
            % 
            % S_21=V_3_base*conj(-I_12);
            % S_23=-S_21-S_L2;
            % 
            % I_23=conj(S_23/V_3_base);
            % 
            % Z_eq_23=S_losses_23/abs(I_23)^2;
            % 
            % T_load.P(1)=real(S_L1);
            % T_load.Q(1)=imag(S_L1);
            % T_load.P(2)=real(S_L2);
            % T_load.Q(2)=imag(S_L2);
            % T_load.P(3)=real(S_L3);
            % T_load.Q(3)=imag(S_L3);
            % 
            % T_load(1,:)=[];
            % 
            % T_NET.R(1)=real(Z_eq_12);
            % T_NET.X(1)=imag(Z_eq_12);
            % 
            % T_NET.R(2)=real(Z_eq_23);
            % T_NET.X(2)=imag(Z_eq_23);

            alpha_p=var_grid(1);
            alpha_q=var_grid(2);
            
            beta_2p=var_grid(3);
            beta_2q=var_grid(4);

            Z_eq_12=(1-info.network.fact_line)*info.network.R_12+(1-info.network.fact_line)*info.network.X_12*1i;
            S_12=info.network.S_net_base;
            I_12=conj(S_12/info.network.V_2_base);

            S_iny_2=info.network.S_iny_base-abs(I_12)^2*Z_eq_12;

            V_3_base=-I_12*Z_eq_12+info.network.V_2_base;

            S_load=alpha_p*real(S_iny_2)+alpha_q*imag(S_iny_2)*1i;
            S_losses=S_iny_2-S_load;

            S_L2=beta_2p*real(S_load)+beta_2q*imag(S_load)*1i;
            S_L3=S_load-S_L2;

            S_21=V_3_base*conj(-I_12);
            S_23=-S_21-S_L2;

            I_23=conj(S_23/V_3_base);

            Z_eq_23=S_losses/abs(I_23)^2;
            % 
            % T_load.P(1)=0;
            % T_load.Q(1)=0;
            T_load.P(1)=real(S_L2);
            T_load.Q(1)=imag(S_L2);
            T_load.P(2)=real(S_L3);
            T_load.Q(2)=imag(S_L3);

            % T_load(1,:)=[];
            T_NET.R(1)=real(Z_eq_12);
            T_NET.X(1)=imag(Z_eq_12);

            T_NET.R(2)=real(Z_eq_23);
            T_NET.X(2)=imag(Z_eq_23);



    end

    T_VSC.mode=info.network.VSC_mode;

    writetable(T_NET, info.newmainfile,'Sheet','AC-NET')
    writetable(T_load, info.newmainfile,'Sheet','load')
    writetable(T_VSC, info.newmainfile,'Sheet','VSC')
    
    if bool_reload_vsc
        T_VSC_data = readtable(info.newvscfile,'Sheet',info.network.VSC_mode{:});
        T_VSC_data{1,3:end} = T_VSC_data{1,3:end}.*10.^(var_vsc);
        writetable(T_VSC_data, info.newvscfile,'Sheet',info.network.VSC_mode{:});
    end 
    


end

