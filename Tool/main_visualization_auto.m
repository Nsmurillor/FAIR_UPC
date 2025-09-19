close all;
clearvars
clc
path(pathdef)
addpath(genpath(pwd))

fact_line=0.01;
R_12=0.000575259116172393;
X_12=0.000293244885684409;
P_sel=0.01;
Q_sel=0;

fig_folder='05_figures\'; 
res_folder='02_results\';

excel_plt    = ['plot_excel.xlsx']; 
sheets_plt = sheetnames(excel_plt);


T_PLT     = readtable(excel_plt,'Sheet','main_plot');
T_PLT_run = T_PLT(strcmp(T_PLT.run_plt,'true'),:);


for ii_run = 1:height(T_PLT_run)
    close all
    T_PLT_ii = T_PLT_run(ii_run,:);

    if strcmp(T_PLT_ii.save_data,'true')
        folder_results=[fig_folder,T_PLT_ii.plot_case{:},'\'];
    else
        folder_results=[fig_folder,'temp\'];
    end

    if exist(folder_results)
        rmdir(folder_results,'s')
    end

    mkdir(folder_results)
    
    bool_sheets=sum(ismember(sheets_plt,T_PLT_ii{:,{'cases_sheet','config_sheet','style_sheet'}})) == 3;
    if bool_sheets
        bool_aux=any(ismember(T_PLT_ii{:,{'PF','DSI','SVD','Y_bode'}},'true'));
        bool_EMT=strcmp(T_PLT_ii.EMT,'true');
        bool_DSI=strcmp(T_PLT_ii.DSI,'true');

        bool_zoom=strcmp(T_PLT_ii.zoom_EMT,'true');

        T_Cases  = readtable(excel_plt,'Sheet',T_PLT_ii.cases_sheet{:});
        T_Style  = readtable(excel_plt,'Sheet',T_PLT_ii.style_sheet{:});
        T_config = readtable(excel_plt,'Sheet',T_PLT_ii.config_sheet{:});

        T_Cases.File_names_aux=strcat(T_Cases.converter,'_aux_',T_Cases.main_name,'.mat');
        T_Cases.File_names_EMT=strcat(T_Cases.converter,'_EMT_',T_Cases.main_name,'.mat');

        T_Cases_run = T_Cases(strcmp(T_Cases.plot,'true'),:);
        if bool_EMT
            T_Cases_run.e_vq1 = nan(height(T_Cases_run),1);
            T_Cases_run.e_vd1 = nan(height(T_Cases_run),1);
            T_Cases_run.e_iq_th = nan(height(T_Cases_run),1);
            T_Cases_run.e_id_th = nan(height(T_Cases_run),1);
        end

        if bool_DSI
            T_Cases_run.DSI_H2_norm_tr = nan(height(T_Cases_run),1);
            T_Cases_run.DSI_H2_norm_SVD = nan(height(T_Cases_run),1);
           
        end


        clear info

        info.path.results=res_folder;

        info.network.fact_line = fact_line;
        info.network.R_12 = R_12;
        info.network.X_12 = X_12;
        info.network.P_sel = P_sel;
        info.network.Q_sel = Q_sel;
        

        for jj_case=1:height(T_Cases_run)
            T_Cases_jj=T_Cases_run(jj_case,:);
            bool_EMT_name=isfile(['02_results/',T_Cases_jj.File_names_EMT{:}]);
            bool_aux_name=isfile(['02_results/',T_Cases_jj.File_names_aux{:}]);

            bool_plot = strcmp(T_Cases_jj.plot,'true');
            bool_force = strcmp(T_Cases_jj.force_sim,'true');
            
            if bool_aux & bool_plot & (bool_force | ~bool_aux_name)  
                info.T_Case=T_Cases_jj;
                info.case='AUX';
                run sim_case
            end

            if bool_EMT & bool_plot & (bool_force | ~bool_EMT_name)  
                info.T_Case=T_Cases_jj;
                info.case='EMT';
                run sim_case
            end

        end
        
        clear info 

        all_plots={'PF','DSI','SVD','Y_bode','EMT'};        
        
        base_cases=unique(T_Cases.idx_base);
        all_cases=unique(T_Cases.idx);

        % Load_files

        info.AUX=cell(1,max(all_cases));
        info.EMT=cell(1,max(all_cases));

        for jj_case=1:height(T_Cases_run) 
            T_Cases_jj = T_Cases_run(jj_case,:);

            if bool_EMT
                info.EMT{T_Cases_jj.idx}=load(['02_results/',T_Cases_jj.File_names_EMT{:}]); 
            end            

            if bool_aux
                info.AUX{T_Cases_jj.idx}=load(['02_results/',T_Cases_jj.File_names_aux{:}]); 
            end
        end


        if bool_aux

            w_data=2*pi*logspace(log10(T_PLT_ii.f_min),log10(T_PLT_ii.f_max),T_PLT_ii.size);
            w_marker=2*pi*logspace(log10(T_PLT_ii.f_min),log10(T_PLT_ii.f_max),30);
    
            idx_mkr = zeros(size(w_marker));
    
            for i = 1:length(w_marker)
                [~, idx] = min(abs(w_data - w_marker(i)));
                idx_mkr(i) = idx;
            end
            
            idx_mkr(end)=[];

            idx_mkr_del=floor([0:height(T_Cases_run)-1].*(min(diff(idx_mkr)))/height(T_Cases_run));

        end


        for jj=1:size(all_plots,2)            
            run init_set_figures.m
        end



        for jj_case=1:height(T_Cases_run)

            T_Cases_jj = T_Cases_run(jj_case,:);
            T_Style_jj =  T_Style(T_Style.idx_style==T_Cases_jj.idx_style,:);

            if isempty(str2num(T_Style_jj.Color{:}))
                color_jj=T_Style_jj.Color{:};
            else
                color_jj=str2num(T_Style_jj.Color{:});
            end     

            if bool_EMT; EMT=info.EMT{T_Cases_jj.idx}; end            
            
            if bool_aux; AUX=info.AUX{T_Cases_jj.idx}; end

            run plot_figures.m

        end

        run final_set_figures.m

    else
        disp(['Some auxiliary sheets not found: ',T_PLT_ii.plot_case{:}])
    end
   
end

%%

close all