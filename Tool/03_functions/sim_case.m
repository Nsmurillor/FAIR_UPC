% function sim_case(info)

    caseName=info.T_Case.main_name{:};
    fanals = 2; 
    shared_power = 0;
    run set_file_names.m
    run read_data.m
    
    
    T_VSC.mode(:)=info.T_Case.converter;
    
    
    run preprocess_data.m
    run get_parameters.m
    
    set_breaker_state('line',1,'close')
    
    
    if true
    
    num_nodes=length(unique([T_NET.bus_from,T_NET.bus_to]));
    
    if    num_nodes>2
        T_NET(1,:).B=0;
        T_load.bus=T_load.bus+1;
        T_VSC.bus=T_VSC.bus+1;        
        
        T_NET=[T_NET(1,:);T_NET];
        T_NET.bus_from(2:end)=T_NET.bus_from(2:end)+1;
        T_NET.bus_to(2:end)=T_NET.bus_to(2:end)+1;
        T_NET.number=[1:height(T_NET)]';
        T_NET.R(1)=info.network.fact_line*info.network.R_12;
        T_NET.X(1)=info.network.fact_line*info.network.X_12;
        T_NET.R(2)=(1-info.network.fact_line)*info.network.R_12;
        T_NET.X(2)=(1-info.network.fact_line)*info.network.X_12;
        
        T_PF=[T_PF;T_PF(end,:)];
        T_PF.bus(end)=T_PF.bus(end)+1;
        
        T_Aggr=[T_Aggr(1,:);T_Aggr];
        T_Aggr(2:end,:)=T_Aggr(2:end,:)+1;
    
    elseif num_nodes==2
    
        T_NET(1,:).B=0;
        T_load.bus=T_load.bus+1;
        T_VSC.bus=T_VSC.bus+1;     
        
        T_NET=[T_NET(1,:);T_NET];
        T_NET.bus_from(2:end)=T_NET.bus_from(2:end)+1;
        T_NET.bus_to(2:end)=T_NET.bus_to(2:end)+1;
        T_NET.number=[1:height(T_NET)]';
        T_NET.R(1)=info.network.fact_line*info.network.R_12;
        T_NET.X(1)=info.network.fact_line*info.network.X_12;
        
        T_PF=[T_PF;T_PF(end,:)];
        T_PF.bus(end)=T_PF.bus(end)+1;
        
        T_Aggr=[T_Aggr(1,:);T_Aggr];
        T_Aggr(2:end,:)=T_Aggr(2:end,:)+1;
    else
        T_NET(1,:)={1,1,1,1,2,info.network.fact_line*info.network.R_12,info.network.fact_line*X_12,0,1}
        
        T_load.bus=T_load.bus+1;
        T_VSC.bus=T_VSC.bus+1;
        
        T_PF=[T_PF;T_PF(end,:)];
        T_PF.bus(end)=T_PF.bus(end)+1;
        
        T_Aggr=[T_Aggr(1,:);T_Aggr];
        T_Aggr(2:end,:)=T_Aggr(2:end,:)+1;
    
    end
    end
    if true
    T_load=[T_load(1,:);T_load];
    T_load.bus(1)=1;
    T_load.P(1)=info.network.P_sel;
    T_load.Q(1)=info.network.Q_sel;
    T_load.number=[1:height(T_load)]';
    end

    switch info.case
        case 'AUX'
            T_NET(1,:).B=0;
        case 'EMT'
            T_NET(1,:).B=1e-5;
            T_TH.R=0.005; %SCR=20
            T_TH.X=0.05;
    end

    run PF_results.m;
    run update_OP.m
    run delta_slack_acdc.m
    run generate_NET_with_Qneg.m
    run generate_elements.m

    switch info.case
        case 'AUX'
            inputNames = {};

            for idx = 1:length(l_blocks)
                ss_block = l_blocks{idx};
                inputNames = cat(2, inputNames, ss_block.InputName');
            end
            
            inputNames = sort(unique(inputNames)); % display BOTH internal and external outputs
            input = inputNames;
        
            run select_all_outputs.m
        
            ss_sys_all = connect(l_blocks{:}, input, output);
        
            sys_data.T_Aggr=T_Aggr;
            sys_data.T_PF=T_PF;
            sys_data.l_blocks=l_blocks;
            sys_data.system_all=ss_sys_all; 
            save([info.path.results,info.T_Case.File_names_aux{:}],'sys_data')

        case 'EMT'

            input = {'TH1.vnq'};
            run select_all_outputs.m
            ss_sys = connect(l_blocks{:}, input, output);
        
            path_fil_slx=['00_tool\Non Linear Models\models\' nonlinear '.slx'];
        
            tokens = regexp(path_fil_slx, '[^\\\/]*?(?=\.slx$)', 'match');
            
            fileName = tokens{1};
        
        
            if exist(fileName, 'file')==4
                close_system(fileName, 0);
            end

            if isfile(path_fil_slx) 
                delete(path_fil_slx)
            end

            

            newSys=true;
        
            run NET_initialization.m           
            run NET_layout_FORCE.m 
            run dependent_states.m
            run param_nonlinear.m

            tic
            out_nolin = sim(nonlinear);
            time_simnolin=toc;
            MsgBoxH = findall(0,'Type','figure','Name','Initial state conflict');
            close(MsgBoxH);

            if exist(fileName, 'file')==4
                close_system(fileName, 0);
            end
            if isfile(path_fil_slx) 
                delete(path_fil_slx)
            end
            if exist( [nonlinear,'.slxc'], 'file')
                delete([nonlinear,'.slxc']);
            end

            sys_data.T_Aggr=T_Aggr;
            sys_data.T_PF=T_PF;
            sys_data.l_blocks=l_blocks;
            sys_data.system_all=ss_sys; 
            sys_data.simnolin=out_nolin;
            sys_data.simtime.nolin=time_simnolin;

            save([info.path.results,info.T_Case.File_names_EMT{:}],'sys_data')
    end
    

% end


