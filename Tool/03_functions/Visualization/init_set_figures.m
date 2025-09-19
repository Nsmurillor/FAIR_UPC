bool_jj=strcmp(T_PLT_ii.(all_plots{jj}),'true');
info.figs.bool.(all_plots{jj})=bool_jj;

Position_jj=str2num(T_config.Position{strcmp(T_config.plot,all_plots{jj})});
if strcmp(all_plots{jj},'PF') & info.figs.bool.(all_plots{jj})
    
    hdl_all=cell(1,max(all_cases));
    hdl_ax_all=cell(1,max(all_cases));
    
    for kk_idx=1:length(base_cases)
        Base_case_name = T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)};
        hdl_kk=figure('Color',[1,1,1],'Position',Position_jj,...
            'Name',[T_PLT_ii.plot_case{:},' - ',all_plots{jj},' - ',...
            Base_case_name]);
        set(hdl_kk,'defaultAxesTickLabelInterpreter','latex');  
        set(hdl_kk,'defaulttextinterpreter','latex');
        set(hdl_kk,'defaultLegendInterpreter','latex');
        set(hdl_kk, 'defaultLineLineWidth', 1);
        set(hdl_kk, 'defaultAxesFontSize', 14);                
        set(hdl_kk, 'defaultTextFontSize', 14);

        hdl_all{base_cases(kk_idx)}=hdl_kk;

        ax_mag=subplot(2,2,1);
        hold(ax_mag,'on')
        grid(ax_mag,"on")
        set(ax_mag,'LineWidth',1,...
            'Fontsize',13);
        ylabel(ax_mag,' Voltage [pu]'); xlabel(ax_mag,'Nodes');                     

        ax_ph=subplot(2,2,2);
        hold(ax_ph,'on')
        grid(ax_ph,"on")
        set(ax_ph,'LineWidth',1,...
            'Fontsize',13);
        ylabel(ax_ph,' Phase [deg]'); xlabel(ax_ph,'Nodes');

        ax_emag=subplot(2,2,3);
        hold(ax_emag,'on')
        grid(ax_emag,"on")
        set(ax_emag,'LineWidth',1,...
            'Fontsize',13);
        ylabel(ax_emag,' Error [pu]'); xlabel(ax_emag,'Nodes');                     

        ax_eph=subplot(2,2,4);
        hold(ax_eph,'on')
        grid(ax_eph,"on")
        set(ax_eph,'LineWidth',1,...
            'Fontsize',13);
        ylabel(ax_eph,' Error [deg]'); xlabel(ax_eph,'Nodes');
        
        
        hdl_ax_all{base_cases(kk_idx)}={ax_mag,ax_ph,ax_emag,ax_eph};


    end
    info.fhdl.(all_plots{jj})=hdl_all;
    info.Fig.(all_plots{jj}).ax=hdl_ax_all; 
end


if ismember(all_plots{jj},{'DSI','SVD'}) & info.figs.bool.(all_plots{jj})                  
    hdl_jj=figure('Color',[1,1,1],'Position',Position_jj,'Name',[T_PLT_ii.plot_case{:},' - ',all_plots{jj}]);
    
    set(hdl_jj,'defaultAxesTickLabelInterpreter','latex');  
    set(hdl_jj,'defaulttextinterpreter','latex');
    set(hdl_jj,'defaultLegendInterpreter','latex');
    set(hdl_jj, 'defaultLineLineWidth', 1);
    set(hdl_jj, 'defaultAxesFontSize', 14);                
    set(hdl_jj, 'defaultTextFontSize', 14);
    
    info.fhdl.(all_plots{jj})=hdl_jj;
    ax_jj=subplot(1,1,1);
    hold(ax_jj,"on");
    set(ax_jj,...
        'XScale','log', ...
        'YScale','log' ...
    );
    ylabel(ax_jj,all_plots{jj}); 
    xlabel(ax_jj,'Frequency [Hz]');
    grid(ax_jj,'on')
    xlim(ax_jj,[T_PLT_ii.f_min,T_PLT_ii.f_max])
    ax_jj.MinorGridColor = [1,1,1]*0.5;
    info.Fig.(all_plots{jj}).ax=ax_jj;                    
end

if strcmp(all_plots{jj},{'Y_bode'}) & info.figs.bool.(all_plots{jj})
    
    hdl_all=cell(1,4);
    hdl_ax_all=cell(1,4);
    Y_name={'qq','qd','dq','dd'};
    
    for kk_idx=1:4
        hdl_kk=figure('Color',[1,1,1],'Position',Position_jj,'Name',[T_PLT_ii.plot_case{:},' - ',all_plots{jj},'_',Y_name{kk_idx}]);
    
        set(hdl_kk,'defaultAxesTickLabelInterpreter','latex');  
        set(hdl_kk,'defaulttextinterpreter','latex');
        set(hdl_kk,'defaultLegendInterpreter','latex');
        set(hdl_kk, 'defaultLineLineWidth', 1);
        set(hdl_kk, 'defaultAxesFontSize', 14);                
        set(hdl_kk, 'defaultTextFontSize', 14);
        
        hdl_all{kk_idx}=hdl_kk;

        ax_mag=subplot(2,1,1);
        hold(ax_mag,'on')
        grid(ax_mag,"on")
        set(ax_mag,'LineWidth',1,...
            'Fontsize',14,...
            'XScale','log' ...
            );
        ylabel(ax_mag,' Mag [dB]'); xlabel(ax_mag,'Frequency (Hz)');
        ax_mag.MinorGridColor = [1,1,1]*0.5;

        ax_ph=subplot(2,1,2);
        hold(ax_ph,'on')
        grid(ax_ph,"on")
        set(ax_ph,'LineWidth',1,...
            'Fontsize',14,...
            'XScale','log' ...
            );
        ylabel(ax_ph,' Phase [deg]'); xlabel(ax_ph,'Frequency (Hz)');
        ax_ph.MinorGridColor = [1,1,1]*0.5;
        
        hdl_ax_all{kk_idx}={ax_mag,ax_ph};

    end
    info.fhdl.(all_plots{jj})=hdl_all;
    info.Fig.(all_plots{jj}).ax=hdl_ax_all; 
end

if strcmp(all_plots{jj},{'EMT'}) & info.figs.bool.(all_plots{jj})

    hdl_all=cell(1,max(all_cases));
    hdl_ax_all=cell(1,max(all_cases));

    x_zoom = [0.04,0.5];
    
    
    dx=diff(x_zoom);
    
    if bool_zoom
        x2_zoom = [0.45,0.5];
        dx2=diff(x2_zoom);
        
        
    end



    Figure_EMT_info=[{'vq1'},{'$v_{1}^{q}$ [pu]'},{'$\int (\Delta v_{1}^{q})^{2} dt$'},{'vn1q'};...
                     {'vd1'},{'$v_{1}^{d}$ [pu]'},{'$\int (\Delta v_{1}^{d})^{2} dt$'},{'vn1d'};...
                     {'iq_th'},{'$i_{TH}^{q}$ [pu]'},{'$\int (\Delta i_{TH}^{q})^{2} dt$'},{'I_TH1_q'};...
                     {'id_th'},{'$i_{TH}^{d}$ [pu]'},{'$\int (\Delta i_{TH}^{d})^{2} dt$'},{'I_TH1_d'}];

    for kk_idx=1:length(base_cases)
        Base_case_name = T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)};
        
        hdl_all_in=cell(1,height(Figure_EMT_info));
        hdl_all_ax_in=cell(1,height(Figure_EMT_info));
        
        if bool_zoom
            bool_GFOL_base=strcmp(T_Cases.converter(T_Cases.idx==base_cases(kk_idx)),'GFOL');

            if bool_GFOL_base
                yp_zoom=[0.65,0.8,0.65,0.65];
            else
                yp_zoom=[0.65,0.8,0.8,0.65];
            end
        end
            
        
        for hh_idx=1:height(Figure_EMT_info)
            hdl_hh=figure('Color',[1,1,1],'Position',Position_jj,...
                'Name',[T_PLT_ii.plot_case{:},' - EMT ',Figure_EMT_info{hh_idx},' - ',...
                Base_case_name]);
            set(hdl_hh,'defaultAxesTickLabelInterpreter','latex');  
            set(hdl_hh,'defaulttextinterpreter','latex');
            set(hdl_hh,'defaultLegendInterpreter','latex');
            set(hdl_hh, 'defaultLineLineWidth', 1);
            set(hdl_hh, 'defaultAxesFontSize', 14);                
            set(hdl_hh, 'defaultTextFontSize', 14);

            hdl_all_in{hh_idx}=hdl_hh;

            ax_mag=subplot(2,1,1);
            hold(ax_mag,'on')
            grid(ax_mag,"on")
            set(ax_mag,'LineWidth',1,...
                'Fontsize',13);
            ylabel(ax_mag,Figure_EMT_info{hh_idx,2}); xlabel(ax_mag,'Time [s]'); 
            xlim(ax_mag,x_zoom)

            

            if bool_zoom
                ax_mag_zoom=axes(hdl_hh,'Position', [0.3,yp_zoom(hh_idx),0.5,0.1]);
                box(ax_mag_zoom,"on")
                hold(ax_mag_zoom,"on")
                set(ax_mag_zoom, 'XTick', []);
                set(ax_mag_zoom, 'YTick', []);                
            else
                ax_mag_zoom=[];

            end



            ax_e=subplot(2,1,2);
            hold(ax_e,'on')
            grid(ax_e,"on")
            set(ax_e,'LineWidth',1,...
                'Fontsize',13);
            ylabel(ax_e,Figure_EMT_info{hh_idx,3}); xlabel(ax_e,'Time [s]'); 
            xlim(ax_e,x_zoom)
            
            
            hdl_all_ax_in{hh_idx}={ax_mag,ax_e,ax_mag_zoom};
        end

        hdl_all{base_cases(kk_idx)}=hdl_all_in;
        hdl_ax_all{base_cases(kk_idx)}=hdl_all_ax_in;

    end
    info.fhdl.(all_plots{jj})=hdl_all;
    info.Fig.(all_plots{jj}).ax=hdl_ax_all; 

end