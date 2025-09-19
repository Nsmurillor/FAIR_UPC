if info.figs.bool.PF

   for kk_idx=1: length(base_cases)
       
       lgnd=legend(info.Fig.PF.ax{base_cases(kk_idx)}{1},...
           'Location','best','NumColumns',6,'FontSize',12);

       set(info.Fig.PF.ax{base_cases(kk_idx)}{1},...
           'Position',info.Fig.PF.ax{base_cases(kk_idx)}{1}.Position-[0,0.07,0,0])
       set(info.Fig.PF.ax{base_cases(kk_idx)}{2},...
           'Position',info.Fig.PF.ax{base_cases(kk_idx)}{2}.Position-[0,0.07,0,0])
       set(info.Fig.PF.ax{base_cases(kk_idx)}{3},...
           'Position',info.Fig.PF.ax{base_cases(kk_idx)}{3}.Position-[0,0.03,0,0])
       set(info.Fig.PF.ax{base_cases(kk_idx)}{4},...
           'Position',info.Fig.PF.ax{base_cases(kk_idx)}{4}.Position-[0,0.03,0,0])

       set(lgnd,'Position',[0.3103    0.8909    0.4099    0.0624])

       saveas(info.fhdl.PF{base_cases(kk_idx)},...
           [folder_results,'PF_',T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)},'.png']);
       savefig(info.fhdl.PF{base_cases(kk_idx)},...
           [folder_results,'PF_',T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)},'.fig']);
        

   end
end

if info.figs.bool.DSI
    legend(info.Fig.DSI.ax,'Location','best')
    saveas(info.fhdl.DSI, [folder_results,'DSI.png']);
    savefig(info.fhdl.DSI,[folder_results,'DSI.fig']);
    if size(base_cases,1)==1
        H2_norm_results=T_Cases_run(:,[{'idx'},{'main_name'},{'DSI_H2_norm_tr'},{'DSI_H2_norm_SVD'}]);
        writetable(H2_norm_results, [folder_results,'Summary','_',T_Cases.main_name{T_Cases.idx==base_cases},'.xlsx'], 'Sheet', 'H2_norm');
    end        
end
if info.figs.bool.SVD
    legend(info.Fig.SVD.ax,'Location','best')
    saveas(info.fhdl.SVD, [folder_results,'SVD.png']);
    savefig(info.fhdl.SVD,[folder_results,'SVD.fig']);
end

if info.figs.bool.Y_bode
    Y_name={'qq','qd','dq','dd'};
    for kk_idx=1:4
        lgnd=legend(info.Fig.Y_bode.ax{kk_idx}{1},'Location','best','NumColumns',4,'FontSize',12);
        set(info.Fig.Y_bode.ax{kk_idx}{1},'Position',[0.1300 0.5838-0.07 0.7750 0.3412])
        set(lgnd,'Position',[0.0674 0.9395-0.04 0.9062 0.0730])
        
        saveas(info.fhdl.Y_bode{kk_idx}, [folder_results,'bode_',Y_name{kk_idx},'.png']);
        savefig(info.fhdl.Y_bode{kk_idx},[folder_results,'bode_',Y_name{kk_idx},'.fig']);
    end

end

if info.figs.bool.EMT
   for kk_idx=1: length(base_cases)

       for hh_idx=1:height(Figure_EMT_info)
           
          lgnd=legend(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1},...
               'Location','best','NumColumns',3,'FontSize',12);
    
           set(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1},...
               'Position',info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.Position-[0,0.07,0,0])
           set(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{2},...
               'Position',info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{2}.Position-[0,0.03,0,0])

           set(lgnd,'Position',[0.2    0.8909    0.6    0.08])

           if bool_zoom
               set(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3},...
                   'Position',info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.Position-[0,0.07,0,0])

                x_lim_1=info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.XLim;
                y_lim_1=info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.YLim;
                dely_1=y_lim_1(2)-y_lim_1(1);
                y_lim_1=y_lim_1+1*[-dely_1,dely_1];
                info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.YLim=info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.YLim+0.5*[-dely_1,dely_1];
                Points_1_rect=[x_lim_1(1) y_lim_1(1); x_lim_1(2)-x_lim_1(1) y_lim_1(2)-y_lim_1(1)];
        
                rectangle(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1},'Position',[Points_1_rect(1,:),Points_1_rect(2,:)],...
                    'LineStyle','--')
        
        
                coor_1_rec=[Points_1_rect(1,:);Points_1_rect(1,:)+Points_1_rect(2,:)];
                Pos_1_box=[info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.InnerPosition(1:2);...
                    info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.InnerPosition(3:4)+info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{3}.InnerPosition(1:2)];
                
                coor_1_rec=[coor_1_rec;coor_1_rec(2,1),coor_1_rec(1,2);coor_1_rec(1,1),coor_1_rec(2,2)];
                
                Pos_1_box=[Pos_1_box;Pos_1_box(2,1),Pos_1_box(1,2);Pos_1_box(1,1),Pos_1_box(2,2)];
                
                coor_1_grid=[info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.XLim(1),info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.YLim(1);info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.XLim(2),info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.YLim(2)]; 
                Pos_1_grid=[info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.InnerPosition(1:2);info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.InnerPosition(3:4)+info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1}.InnerPosition(1:2)];
                
                delx1_pos_grid=Pos_1_grid(2,1)-Pos_1_grid(1,1);
                delx1_coor=coor_1_grid(2,1)-coor_1_grid(1,1);
                
                dely1_pos_grid=Pos_1_grid(2,2)-Pos_1_grid(1,2);
                dely1_coor=coor_1_grid(2,2)-coor_1_grid(1,2);
                
        
                for kk=1:4
                    pos_box_i_kk=Pos_1_box(kk,:);
                    pos_rec_kk=coor_1_rec(kk,:);
        
                    coor_x_1=coor_1_grid(1,1)+(pos_box_i_kk(1)-Pos_1_grid(1,1))*(delx1_coor)/(delx1_pos_grid);
                    coor_y_1=coor_1_grid(1,2)+(pos_box_i_kk(2)-Pos_1_grid(1,2))*(dely1_coor)/(dely1_pos_grid);
        
                    plot(info.Fig.EMT.ax{base_cases(kk_idx)}{hh_idx}{1},[coor_x_1,pos_rec_kk(1)],[coor_y_1,pos_rec_kk(2)],'Color','k',...
                        'LineStyle',':','HandleVisibility','off')
                end
           end

    
           saveas(info.fhdl.EMT{base_cases(kk_idx)}{hh_idx},...
               [folder_results,'EMT_',Figure_EMT_info{hh_idx,1},'_',T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)},'.png']);
           saveas(info.fhdl.EMT{base_cases(kk_idx)}{hh_idx},...
               [folder_results,'EMT_',Figure_EMT_info{hh_idx,1},'_',T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)},'.fig']);
       end

        EMT_results=T_Cases_run(:,[{'idx'},{'main_name'},{'e_vq1'},{'e_vd1'},{'e_iq_th'},{'e_id_th'}]);
        writetable(EMT_results, [folder_results,'Summary','_',T_Cases.main_name{T_Cases.idx==base_cases(kk_idx)},'.xlsx'], 'Sheet', 'EMT_errors');

   end

end