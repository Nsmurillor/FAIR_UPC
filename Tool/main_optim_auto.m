close all;
clearvars
clc
path(pathdef)
addpath(genpath(pwd))

% set(groot,'defaultAxesTickLabelInterpreter','latex');  
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% Parameters special modificaction

fact_line=0.01;
R_12=0.000575259116172393;
X_12=0.000293244885684409;
P_sel=0.01;
Q_sel=0;

opt_folder='04_results_opt\'; 

% Read Excel Data iterations

excel_opt    = ['opt_excel_01.xlsx']; 
sheets_opt = sheetnames(excel_opt);

T_OPT     = readtable(excel_opt,'Sheet','main_optimization');
T_OPT_run = T_OPT(strcmp(T_OPT.run_case,'true'),:);


Times_xd=[];
for ii_run = 1:height(T_OPT_run)
    tic
    
    T_OPT_ii = T_OPT_run(ii_run,:);

    if strcmp(T_OPT_ii.save_data,'true')
        folder_results=['04_results_opt\',T_OPT_ii.test_name{:},'\'];
    else
        folder_results='04_results_opt\temp\';
    end

    if exist(folder_results)
        rmdir(folder_results,'s')
    end

    mkdir(folder_results)
    

    excel_base = [T_OPT_ii.main_file,{'System_ref'}];
    kron_base = [T_OPT_ii.kron_file,{'System_kron'}];
    excel_1gen_base = [T_OPT_ii.base_file_1gen,{'System_base'}];

    bool_GFOR=strcmp(T_OPT_ii.converter,'GFOR');

    cases_vec = [excel_base;kron_base];

    for ii_case=1:height(cases_vec)
        


        caseName = cases_vec{ii_case,1};
        fanals = 2; 

        shared_power = 0;
        run set_file_names.m
        run read_data.m

        if bool_GFOR 
            T_VSC.mode(:)={'GFOR'};
        end

        run preprocess_data.m
        run get_parameters.m

        set_breaker_state('line',1,'close')
        
        % Add_extra node and load
        T_NET(1,:).B=0;
        T_load.bus=T_load.bus+1;
        T_VSC.bus=T_VSC.bus+1;
        
        T_NET=[T_NET(1,:);T_NET];
        T_NET.bus_from(2:end)=T_NET.bus_from(2:end)+1;
        T_NET.bus_to(2:end)=T_NET.bus_to(2:end)+1;
        T_NET.number=[1:height(T_NET)]';
        T_NET.R(1)=fact_line*R_12;
        T_NET.X(1)=fact_line*X_12;
        T_NET.R(2)=(1-fact_line)*R_12;
        T_NET.X(2)=(1-fact_line)*X_12;
        
        T_PF=[T_PF;T_PF(end,:)];
        T_PF.bus(end)=T_PF.bus(end)+1;
        
        T_Aggr=[T_Aggr(1,:);T_Aggr];
        T_Aggr(2:end,:)=T_Aggr(2:end,:)+1;

        T_load=[T_load(1,:);T_load];
        T_load.bus(1)=1;
        T_load.P(1)=P_sel;
        T_load.Q(1)=Q_sel;
        T_load.number=[1:height(T_load)]';

        run PF_results.m
        run update_OP.m
        run delta_slack_acdc.m
        run generate_NET_with_Qneg.m
        run generate_elements.m


        input={"NET.vn1q","NET.vn1d"};
        output={"NET.iq_1_2","NET.id_1_2"};
        
        l_block_base = l_blocks(:);
        l_block_base([2]) = []; %Load disconection


        switch cases_vec{ii_case,2}

            case 'System_ref'

                [Ybus, Yf, Yt] = makeYbus(mpc);
                net_base=network(T_NET,T_load,T_VSC,T_TH,T_SG,full(Ybus),T_PF.Vm,T_PF.theta);
                net_base=net_base.runcarpintery(); 
                
                S_net_base=-net_base.results.S_tt(2,1);
                S_gen_base=sum(T_VSC.P)+sum(T_VSC.Q)*1i;
                S_iny_base=S_net_base+S_gen_base;
                V_2_base=net_base.results.V(2);

                System_ref=(connect(l_block_base{:},input,output));

            case 'System_kron'

                System_kron=(connect(l_block_base{:},input,output));

        end

    end


    diff_kron=System_ref-System_kron;

    % Optimization setup
    clear info
    info.caseName=excel_1gen_base{:,1};
    
    info.network.V_2_base = V_2_base;
    info.network.S_net_base = S_net_base;
    info.network.S_gen_base =S_gen_base ;
    info.network.S_iny_base= S_iny_base;

    info.network.bool_GFOR = bool_GFOR;
    info.network.fact_line = fact_line;
    info.network.VSC_mode = T_OPT_ii.converter;
    info.network.R_12 = R_12;
    info.network.X_12 = X_12;
    info.network.P_sel = P_sel;
    info.network.Q_sel = Q_sel;
    

    info.System.System_ref=System_ref;
    info.System.System_kron=System_kron;

    info.type=T_OPT_ii.type{:};
    info.version=T_OPT_ii.version;

    info.Data_preload=preload_data(excel_1gen_base{:,1},info);

    bool_limits_sheet=ismember(T_OPT_ii.limits,sheets_opt);

   
   
    if bool_limits_sheet
        info.lim.bool=true;
        lims=readtable(excel_opt,'Sheet',T_OPT_ii.limits{:});
        info.lim.number=height(lims);

        info.lim.w_opt=[];
        info.lim.DSI_kron=[];

        for ii_lim_idx=1:info.lim.number
            ii_lim=lims(ii_lim_idx,:);
            w_opt_ii=2*pi*logspace(log10(ii_lim.lim_inf),log10(ii_lim.lim_sup),ii_lim.size);
            DSI_kron_ii=sigma(diff_kron,w_opt_ii);
            info.lim.w_opt=[info.lim.w_opt;{w_opt_ii}];
            info.lim.DSI_kron=[info.lim.DSI_kron;{DSI_kron_ii(1,:)}];
        end

        info.lim.w_opt_log=cellfun(@(x) log10(x),info.lim.w_opt, 'UniformOutput', false);
        info.lim.DSI_kron_log=cellfun(@(x) log10(x),info.lim.DSI_kron, 'UniformOutput', false);
        info.lim.factor=lims.factor;
        info.lim.bias=lims.bias;
        
    else
        info.lim.number=0;
        info.lim.bool=false;
    end

    bool_line_sheet=ismember(T_OPT_ii.lines,sheets_opt);

    if bool_line_sheet
        info.line.bool=true;

        lines=readtable(excel_opt,'Sheet',T_OPT_ii.lines{:});
        info.line.number=height(lines);

        info.line.w_opt=[];
        info.line.DSI_line=[];
        info.line.DSI_kron=[];


        for ii_line_idx=1:info.line.number
            ii_line=lines(ii_line_idx,:);
            log10_inf=log10(ii_line.f_inf*2*pi);
            log10_sup=log10(ii_line.f_sup*2*pi);
            
            w_opt_ii=logspace(log10_inf,log10_sup,ii_line.size);
            DSI_line_ii=10.^((ii_line.DSI_sup-ii_line.DSI_inf)/(log10_sup-log10_inf).*(log10(w_opt_ii)-log10_inf)+ii_line.DSI_inf);
            DSI_kron_ii=sigma(diff_kron,w_opt_ii);

            info.line.w_opt=[info.line.w_opt;{w_opt_ii}];
            info.line.DSI_line=[info.line.DSI_line;{DSI_line_ii}];
            info.line.DSI_kron=[info.line.DSI_kron;{DSI_kron_ii(1,:)}];
        end

        info.line.w_opt_log=cellfun(@(x) log10(x),info.line.w_opt, 'UniformOutput', false);
        info.line.DSI_line_log=cellfun(@(x) log10(x),info.line.DSI_line, 'UniformOutput', false);
        info.line.DSI_kron_log=cellfun(@(x) log10(x),info.line.DSI_kron, 'UniformOutput', false);
        info.line.factor=lines.factor;
    else
        info.line.number=0;
        info.line.bool=false;
    end

    bool_points_sheet=ismember(T_OPT_ii.points,sheets_opt);

    if bool_points_sheet
        info.point.bool=true;

        points=readtable(excel_opt,'Sheet',T_OPT_ii.points{:});
        info.point.number=height(points);

        info.point.w_opt=[];
        info.point.DSI_kron=[];

        for ii_point_idx=1:info.point.number
            ii_point=points(ii_point_idx,:);
            w_opt_ii=2*pi*ii_point.value;
            DSI_kron_ii=sigma(diff_kron,w_opt_ii);
            info.point.w_opt=[info.point.w_opt;w_opt_ii];
            info.point.DSI_kron=[info.point.DSI_kron;DSI_kron_ii(1,:)];
        end

        info.point.w_opt_log=log10(info.point.w_opt);
        info.point.DSI_kron_log=log10(info.point.DSI_kron);
        info.point.factor=points.factor;
        info.point.bias=points.bias;

    else
        info.point.number=0;
        info.point.bool=false;
    end
    



     net_var=4;


    switch info.type

        case 'Grid'
            nvars = net_var;
            lb = [0.001*ones(1,net_var)];
            ub = [0.999*ones(1,net_var)];

        case 'Electric'
            nvars = 6+net_var;
            lb = [0.001*ones(1,net_var),-0.2*ones(1,6)];
            ub = [0.999*ones(1,net_var),0.2*ones(1,6)];
            
        case 'All'
            nvars = 15+net_var;
            lb = [0.001*ones(1,net_var),-0.2*ones(1,15)];
            ub = [0.999*ones(1,net_var),0.2*ones(1,15)];
    end



    info.plot.bool=true;

    if info.plot.bool
        run set_fig_opt.m
    end

    info.optmode='fmin';
    info.plot.time_init=datetime('now');

    % x_try=[0.1,0.2,0.2,0.2];
    % DSI_optim(x_try,info) 

    objFun = @(x) DSI_optim(x,info); 

    local_opts = optimoptions('fmincon', ...
        'Display','iter', ...
        'PlotFcn','optimplotfval', ...
        'UseParallel', true, ...
        'MaxIterations', T_OPT_ii.fmin_iter, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10);
    
    options = optimoptions('particleswarm', ...
        'Display', 'iter', ...             % Muestra progreso en consola
        'SwarmSize', T_OPT_ii.swarm_size, ...  % Número de partículas
        'PlotFcn',@pswplotbestf,...
        'InertiaRange', [0.3 0.7], ...           % Menor inercia = más explotación
        'SelfAdjustmentWeight', 1.5, ...         % Peso hacia posición propia
        'SocialAdjustmentWeight', 1.0, ...       % Peso hacia el mejor global
        'HybridFcn', {@fmincon,local_opts},...
        'MaxIterations',T_OPT_ii.swarm_iter,... % Número máx de iteraciones
        'UseParallel', true);                   % Parallel computing
    


    [x_opt,fval,exitflag,output_opt,points] = particleswarm(objFun, nvars, lb, ub, options);
    
 
    figs = findall(0, 'Type', 'figure');
    opt_fig = findobj(figs, 'Name', 'Optimization Plot Function');
     

    Base_Path = which([info.caseName,'.xlsx']);
    New_path=fullfile(folder_results, [T_OPT_ii.test_name{:},'.xlsx']);
    copyfile(Base_Path, New_path)

    Base_Path_vsc = which([info.caseName,'_data_vsc.xlsx']);
    New_path_vsc=fullfile(folder_results, [T_OPT_ii.test_name{:},'_data_vsc.xlsx']);
    copyfile(Base_Path_vsc, New_path_vsc)

    info.newmainfile=New_path;
    info.newvscfile=New_path_vsc;
    
    DSI_optim_result(x_opt,info)
    
   
    run save_plots
    
    if info.plot.bool     
        run save_animation
        savefig(info.plot.fig_1,[folder_results,'Final_DSI.fig'])
        exportgraphics(info.plot.fig_1,[folder_results,'Final_DSI.png'], 'Resolution', 300)
    end

    close all;
    time_iter=toc;
    
    Times_xd=[Times_xd,time_iter];

end

