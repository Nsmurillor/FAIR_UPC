function  load_data = preload_data(caseName,info)
    
    fanals = 2; 
    shared_power = 0;
    run set_file_names.m
    run read_data.m

    if info.network.bool_GFOR
        T_VSC.mode(:)={'GFOR'};
    else
        T_VSC.mode(:)={'GFOL'};
    end

    run preprocess_data.m
    run preload_parameters_opt
    
    set_breaker_state('line',1,'close')
    
    num_nodes=length(unique([T_NET.bus_from,T_NET.bus_to]));
    
    switch info.version
        case 1
            T_NET(1,:).B=0;
    
            T_load.bus=T_load.bus+1;
            T_VSC.bus=T_VSC.bus+1;
            T_data_VSC.bus=T_data_VSC.bus+1; 
            T_VSC_base.bus=T_VSC_base.bus+1;
            
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
        case 2
            T_NET(1,:).B=0;
            T_load.bus=T_load.bus+1;
            T_VSC.bus=T_VSC.bus+1;
            T_data_VSC.bus=T_data_VSC.bus+1; 
            T_VSC_base.bus=T_VSC_base.bus+1;
            
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

            T_load=[T_load(1,:);T_load];
            T_load.bus(1)=2;
            T_load.P(1)=0;
            T_load.Q(1)=0;
            T_load.number=[1:height(T_load)]';
    end




    T_load=[T_load(1,:);T_load];
    T_load.bus(1)=1;
    T_load.P(1)=info.network.P_sel;
    T_load.Q(1)=info.network.Q_sel;
    T_load.number=[1:height(T_load)]';


    vars = whos;   % lista de variables locales
    load_data = struct;
    for k = 1:numel(vars)
        nombre = vars(k).name;
        load_data.(nombre) = eval(nombre);
    end

end