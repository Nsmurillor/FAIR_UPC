classdef network
    %NETWORK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        info;
        results;
        Tables;
        
    end
    
    methods
        function obj = network(T_NET,T_LOAD,T_VSC,T_TH,T_SG,Ybus,V,theta)
            obj.Tables.T_VSC=T_VSC;
            obj.Tables.T_TH=T_TH;
            obj.Tables.T_SG=T_SG;     
            obj.Tables.T_LOAD=T_LOAD; 
            obj.Tables.T_NET=T_NET;   
            obj.Tables.T_symp=table();
            obj.info.buses_gen=unique([obj.Tables.T_VSC.bus;obj.Tables.T_TH.bus;obj.Tables.T_SG.bus]);
            
   
   
            obj.info.Y_bus.mpc = Ybus;
            obj.info.nodes = length(Ybus);
            
            obj.results.V=V.*exp(1i*theta*pi/180);
            % obj.Tables.T_symp.old=[1:11,13:34]';
            % obj.Tables.T_symp.old=[1,18,22,25,33]';
            % obj.Tables.T_symp.old=[1,2,18,22,25,33]';
            % obj.Tables.T_symp.old=[1,2,3,4,7,19,23,26,34]';
            % obj.Tables.T_symp.old=[1,2,3,6,17,18,22,25,32,33]';
            % obj.Tables.T_symp.old=[1,2,3,4,6,7,8]';
            % obj.Tables.T_symp.old=[1,2,3,5,6,7]';
            % obj.Tables.T_symp.old=[1,2,4,5]';
            % obj.Tables.T_symp.old=[1,2,4]';
            obj.Tables.T_symp.old=[1,2,3,6,18,22,25,33]';
            % 
            obj.info.newnodes = length(obj.Tables.T_symp.old);
            obj.Tables.T_symp.idx=[1:obj.info.newnodes]';
            
          
            obj.Tables.T_symp.new=[1:obj.info.newnodes]';

            Y1=zeros(obj.info.nodes);
            Y2=zeros(obj.info.nodes);

            variable_names_loads= [{'number'}, {'double'};...
                        {'Area'},{'double'}; ...
                        {'SyncArea'},{'double'};...
                        {'bus'},{'double'};...
                        {'R'},{'double'};...
                        {'X'},{'double'};...
                        {'P'},{'double'};...                            
                        {'Q'},{'double'};...
                        {'type'},{'char'};...
                        {'state'},{'double'}];
           
            obj.Tables.Symp.T_load=table('Size',[0,size(variable_names_loads,1)],... 
	                    'VariableNames', variable_names_loads(:,1),...
	                    'VariableTypes', variable_names_loads(:,2));
            idx_load=1;

            for ii=1:height(T_LOAD) % Add load R to G
                load=T_LOAD(ii,:);
                idx_node=load.bus;
                if strcmp(load.type,'RX')
                    Y1(idx_node,idx_node)=1/(load.R+1i*load.X);
                    Y2(idx_node,idx_node)=1/(load.R+1i*load.X);
                elseif strcmp(load.type,'PQ')
                    % ~ismember(idx_node,obj.info.buses_gen)
                    if true
                        V_ii=abs(obj.results.V(idx_node));
                        Y2(idx_node,idx_node)=(load.P)/V_ii^2-1i*(load.Q)/V_ii^2;
                    else
                        obj.Tables.Symp.T_load.number(idx_load)=idx_load;
                        obj.Tables.Symp.T_load.Area(idx_load)=1;
                        obj.Tables.Symp.T_load.SyncArea(idx_load)=1;
                        obj.Tables.Symp.T_load.bus(idx_load)=obj.Tables.T_symp.new(obj.Tables.T_symp.old==idx_node);
                        obj.Tables.Symp.T_load.R(idx_load)=0;
                        obj.Tables.Symp.T_load.X(idx_load)=0;
                        obj.Tables.Symp.T_load.P(idx_load)=load.P;
                        obj.Tables.Symp.T_load.Q(idx_load)=load.Q;
                        obj.Tables.Symp.T_load.type(idx_load)={'PQ'};
                        obj.Tables.Symp.T_load.state(idx_load)=1;                    
                        idx_load=idx_load+1;

                    end
                end
                
            end
            
            for ii=1:height(T_NET)
                branch=T_NET(ii,:);
                if branch.state==1
                    idx_from=branch.bus_from;
                    idx_to=branch.bus_to;
                    y_ii=1/(branch.R+1i*branch.X);
                    b_ii=1i*branch.B/2;
                    Y1(idx_from,idx_to)=-y_ii;
                    Y1(idx_to,idx_from)=-y_ii;
                    Y1(idx_from,idx_from)=Y1(idx_from,idx_from)+y_ii+b_ii;
                    Y1(idx_to,idx_to)=Y1(idx_to,idx_to)+y_ii+b_ii;

                    Y2(idx_from,idx_to)=-y_ii;
                    Y2(idx_to,idx_from)=-y_ii;
                    Y2(idx_from,idx_from)=Y2(idx_from,idx_from)+y_ii+b_ii;
                    Y2(idx_to,idx_to)=Y2(idx_to,idx_to)+y_ii+b_ii;


                    
                end

            end
            obj.info.Y_bus.Y1=Y1;
            obj.info.Y_bus.Y2=Y2;

        end
        
        function obj=runcarpintery(obj)

            clear Pvec Qvec
            for ii_node=1:obj.info.nodes
                P=0;
                Q=0;
                idx_Y=find(obj.info.Y_bus.Y2(ii_node,:)~=0);
                for jj_node=idx_Y
                    
                    y_iijj=obj.info.Y_bus.Y2(ii_node,jj_node);
                    v_ii=obj.results.V(ii_node);
                    v_jj=obj.results.V(jj_node);
                    P=P+abs(v_ii)*abs(v_jj)*abs(y_iijj)*cos(angle(v_ii)-angle(v_jj)-angle(y_iijj));
                    Q=Q+abs(v_ii)*abs(v_jj)*abs(y_iijj)*sin(angle(v_ii)-angle(v_jj)-angle(y_iijj));                    
                end
                
                
                Pvec(ii_node)=P;
                Qvec(ii_node)=Q;

            end
            obj.results.P=Pvec';
            obj.results.Q=Qvec';
            obj.results.In=obj.info.Y_bus.Y2*obj.results.V;


            I_tt=zeros(obj.info.nodes);
            S_tt=zeros(obj.info.nodes);
            for ii=1:obj.info.nodes
                for jj=1:obj.info.nodes
                    if ii==jj
                        I_tt(ii,jj)=sum(obj.info.Y_bus.Y1(ii,:))*obj.results.V(ii);
                        S_tt(ii,jj)=obj.results.V(ii)*conj(I_tt(ii,jj));
                    else
                        I_tt(ii,jj)=-obj.info.Y_bus.Y1(ii,jj)*(obj.results.V(ii)-obj.results.V(jj));
                        S_tt(ii,jj)=obj.results.V(ii)*conj(I_tt(ii,jj));
                    end
                end
            end

            obj.results.I_tt=I_tt;
            obj.results.S_tt=S_tt;

      
        end
        function obj=simplification(obj)

            V_per=obj.results.V;
            Y_per=obj.info.Y_bus.Y2;
            I_per=obj.results.In;

        

            for ii=1:obj.info.newnodes
                ii_stay=obj.Tables.T_symp.old(ii);
                M_change=eye(obj.info.nodes);
                M_change(:,[ii_stay,ii]) = fliplr(M_change(:,[ii_stay,ii]));
                V_per=M_change*V_per;
                Y_per=M_change*Y_per*M_change;
                I_per=M_change*I_per;                
            end
            Y11=Y_per(1:obj.info.newnodes,1:obj.info.newnodes);
            Y12=Y_per(1:obj.info.newnodes,obj.info.newnodes+1:end);
            Y21=Y_per(obj.info.newnodes+1:end,1:obj.info.newnodes);
            Y22=Y_per(obj.info.newnodes+1:end,obj.info.newnodes+1:end);
            Ynew=Y11-Y12*inv(Y22)*Y21;
            Vnew=V_per(1:obj.info.newnodes);
            Inew=I_per(1:obj.info.newnodes);


            variable_names_ac_net= [{'number'}, {'double'};...
                        {'Area'},{'double'}; ...
                        {'SyncArea'},{'double'};...
                        {'bus_from'},{'double'};...
                        {'bus_to'},{'double'};...
                        {'R'},{'double'};...
                        {'X'},{'double'};...
                        {'B'},{'double'};...
                        {'state'},{'double'}];
            obj.Tables.Symp.T_NET=table('Size',[0,size(variable_names_ac_net,1)],... 
	                    'VariableNames', variable_names_ac_net(:,1),...
	                    'VariableTypes', variable_names_ac_net(:,2));

            idx_branch=1;
            
            for ii=1:obj.info.newnodes

                for jj=ii:obj.info.newnodes
                    if jj>ii
                        if abs(Ynew(ii,jj))>0
                            obj.Tables.Symp.T_NET.number(idx_branch)=idx_branch;
                            obj.Tables.Symp.T_NET.Area(idx_branch)=1;
                            obj.Tables.Symp.T_NET.SyncArea(idx_branch)=1;
                            obj.Tables.Symp.T_NET.bus_from(idx_branch)=ii;
                            obj.Tables.Symp.T_NET.bus_to(idx_branch)=jj;
                            obj.Tables.Symp.T_NET.R(idx_branch)=real(1/-Ynew(ii,jj));
                            obj.Tables.Symp.T_NET.X(idx_branch)=imag(1/-Ynew(ii,jj));
                            obj.Tables.Symp.T_NET.B(idx_branch)=0;
                            obj.Tables.Symp.T_NET.state(idx_branch)=1;                
                            idx_branch=idx_branch+1;
                        end
                    end
    
                end
               
            end
            Y_ground=sum(Ynew);
            idx_load=height(obj.Tables.Symp.T_load)+1;
            
            for ii=1:obj.info.newnodes
                
                    obj.Tables.Symp.T_load.number(idx_load)=idx_load;
                    obj.Tables.Symp.T_load.Area(idx_load)=1;
                    obj.Tables.Symp.T_load.SyncArea(idx_load)=1;
                    obj.Tables.Symp.T_load.bus(idx_load)=ii;
                    % obj.Tables.Symp.T_load.R(idx_load)=real(1/Y_ground(ii));
                    % obj.Tables.Symp.T_load.X(idx_load)=imag(1/Y_ground(ii));                    % 
                    % obj.Tables.Symp.T_load.P(idx_load)=0;
                    % obj.Tables.Symp.T_load.Q(idx_load)=0;
                    obj.Tables.Symp.T_load.R(idx_load)=0;
                    obj.Tables.Symp.T_load.X(idx_load)=0;
                    obj.Tables.Symp.T_load.P(idx_load)=real(abs(Vnew(ii))^2*conj(Y_ground(ii)));
                    obj.Tables.Symp.T_load.Q(idx_load)=imag(abs(Vnew(ii))^2*conj(Y_ground(ii)));
                    obj.Tables.Symp.T_load.type(idx_load)={'PQ'};
                    obj.Tables.Symp.T_load.state(idx_load)=1;
                   
                    idx_load=idx_load+1;
               
            end
            x=12;


        end
        function obj=simplification_electrical(obj)

            V_per=obj.results.V;
            Y_per=obj.info.Y_bus.Y2;
            I_per=obj.results.In;

        

            for ii=1:obj.info.newnodes
                ii_stay=obj.Tables.T_symp.old(ii);
                M_change=eye(obj.info.nodes);
                M_change(:,[ii_stay,ii]) = fliplr(M_change(:,[ii_stay,ii]));
                V_per=M_change*V_per;
                Y_per=M_change*Y_per*M_change;
                I_per=M_change*I_per;                
            end
            Y11=Y_per(1:obj.info.newnodes,1:obj.info.newnodes);
            Y12=Y_per(1:obj.info.newnodes,obj.info.newnodes+1:end);
            Y21=Y_per(obj.info.newnodes+1:end,1:obj.info.newnodes);
            Y22=Y_per(obj.info.newnodes+1:end,obj.info.newnodes+1:end);
            Ynew=Y11-Y12*inv(Y22)*Y21;
            Vnew=V_per(1:obj.info.newnodes);
            Inew=I_per(1:obj.info.newnodes);


            variable_names_ac_net= [{'number'}, {'double'};...
                        {'Area'},{'double'}; ...
                        {'SyncArea'},{'double'};...
                        {'bus_from'},{'double'};...
                        {'bus_to'},{'double'};...
                        {'R'},{'double'};...
                        {'X'},{'double'};...
                        {'B'},{'double'};...
                        {'state'},{'double'}];
            obj.Tables.Symp.T_NET=table('Size',[0,size(variable_names_ac_net,1)],... 
	                    'VariableNames', variable_names_ac_net(:,1),...
	                    'VariableTypes', variable_names_ac_net(:,2));

            idx_branch=1;
            
            for ii=1:obj.info.newnodes

                for jj=ii:obj.info.newnodes
                    if jj>ii
                        if abs(Ynew(ii,jj))>0
                            obj.Tables.Symp.T_NET.number(idx_branch)=idx_branch;
                            obj.Tables.Symp.T_NET.Area(idx_branch)=1;
                            obj.Tables.Symp.T_NET.SyncArea(idx_branch)=1;
                            obj.Tables.Symp.T_NET.bus_from(idx_branch)=ii;
                            obj.Tables.Symp.T_NET.bus_to(idx_branch)=jj;
                            obj.Tables.Symp.T_NET.R(idx_branch)=0;
                            obj.Tables.Symp.T_NET.X(idx_branch)=0;
                            obj.Tables.Symp.T_NET.B(idx_branch)=0;
                            obj.Tables.Symp.T_NET.state(idx_branch)=1;                
                            idx_branch=idx_branch+1;
                        end
                    end
    
                end
               
            end

            Z_iter_={1,[1];...
                     2,[2];...
                     3,[18:21];...
                     4,[3:5];...
                     5,[22:24];...
                     6,[6:17];...
                     7,[25:32]};

            for ii=1:height(obj.Tables.Symp.T_NET)
                idx_ii=Z_iter_{ii,1};
                idx_base=Z_iter_{ii,2};
                R_ii=0;
                X_ii=0;
                for jj=1:length(idx_base)
                    jj_idx=idx_base(jj);
                    R_jj=obj.Tables.T_NET.R(ismember(obj.Tables.T_NET.number,jj_idx));
                    X_jj=obj.Tables.T_NET.X(ismember(obj.Tables.T_NET.number,jj_idx));
                    R_ii=R_jj+R_ii;
                    X_ii=X_jj+X_ii;
                end
                 obj.Tables.Symp.T_NET.R(idx_ii)=R_ii;
                 obj.Tables.Symp.T_NET.X(idx_ii)=X_ii;
            end



            Y_ground=sum(Ynew);
            idx_load=height(obj.Tables.Symp.T_load)+1;
            
            for ii=1:obj.info.newnodes
                
                    obj.Tables.Symp.T_load.number(idx_load)=idx_load;
                    obj.Tables.Symp.T_load.Area(idx_load)=1;
                    obj.Tables.Symp.T_load.SyncArea(idx_load)=1;
                    obj.Tables.Symp.T_load.bus(idx_load)=ii;
                    % obj.Tables.Symp.T_load.R(idx_load)=real(1/Y_ground(ii));
                    % obj.Tables.Symp.T_load.X(idx_load)=imag(1/Y_ground(ii));                    % 
                    % obj.Tables.Symp.T_load.P(idx_load)=0;
                    % obj.Tables.Symp.T_load.Q(idx_load)=0;
                    obj.Tables.Symp.T_load.R(idx_load)=0;
                    obj.Tables.Symp.T_load.X(idx_load)=0;
                    obj.Tables.Symp.T_load.P(idx_load)=0;
                    obj.Tables.Symp.T_load.Q(idx_load)=0;
                    obj.Tables.Symp.T_load.type(idx_load)={'PQ'};
                    obj.Tables.Symp.T_load.state(idx_load)=1;
                   
                    idx_load=idx_load+1;
               
            end

            Y_iter_info={1,1,[];...
                     2,2,[19:21];...
                     3,3,[4:5,23:24];...
                     4,6,[4:5,7:17,26:32];...
                     5,18,[7:17];...
                     6,22,[19:21];...
                     7,25,[23:24];...
                     8,33,[26:32]};
             x=349;
            


            for ii=1:height(obj.Tables.Symp.T_load)
                idx_ii=Y_iter_info{ii,1};
                idx_base=Y_iter_info{ii,2};
                idx_iter=Y_iter_info{ii,3};
                V_ii=abs(obj.results.V(idx_base));
                G_ii=(obj.Tables.T_LOAD.P(ismember(obj.Tables.T_LOAD.bus,idx_base)))/V_ii^2;
                B_ii=-1*(obj.Tables.T_LOAD.Q(ismember(obj.Tables.T_LOAD.bus,idx_base)))/V_ii^2;
                G_sum=0;
                B_sum=0;
                for jj=1:length(idx_iter)
                    jj_idx=idx_iter(jj);
                    V_jj=abs(obj.results.V(jj_idx));
                    G_jj=(obj.Tables.T_LOAD.P(ismember(obj.Tables.T_LOAD.bus,jj_idx)))/V_jj^2;
                    B_jj=-1*(obj.Tables.T_LOAD.Q(ismember(obj.Tables.T_LOAD.bus,jj_idx)))/V_jj^2;

                    G_sum=G_sum+G_jj;
                    B_sum=B_sum+B_jj;
                end
                G_total=G_ii+G_sum/2;
                B_total=B_ii+B_sum/2;

                obj.Tables.Symp.T_load.P(ii)=real(abs(V_ii)^2*conj(G_total+1i*B_total));
                obj.Tables.Symp.T_load.Q(ii)=imag(abs(V_ii)^2*conj(G_total+1i*B_total));
            end
            x=349;
                    % 
                    % obj.Tables.Symp.T_load.P(idx_load)=real(abs(Vnew(ii))^2*conj(Y_ground(ii)));
                    % obj.Tables.Symp.T_load.Q(idx_load)=imag(abs(Vnew(ii))^2*conj(Y_ground(ii)));

        end
        function obj=Gen_aggregation(obj)
            Agg_nodes=[5,6]; % Nodos generadores simplificar
            node_ii=3; % Conexión a dos lineas
            S_load=0;
            S_gen=0;
            S_ii_all=0;
            Losses_all=0;
            for jj=1:length(Agg_nodes)
                node_jj=Agg_nodes(jj);
                P_loadjj=obj.Tables.T_LOAD.P(obj.Tables.T_LOAD.bus==node_jj);
                Q_loadjj=obj.Tables.T_LOAD.Q(obj.Tables.T_LOAD.bus==node_jj);
                P_genjj=obj.Tables.T_VSC.P(obj.Tables.T_VSC.bus==node_jj);
                Q_genjj=obj.Tables.T_VSC.Q(obj.Tables.T_VSC.bus==node_jj);
                S_ii_jj=obj.results.S_tt(node_ii,node_jj);
                S_jj_ii=obj.results.S_tt(node_jj,node_ii);
                Losses_all=Losses_all+S_ii_jj+S_jj_ii;
                S_ii_all=S_ii_all+S_ii_jj;
                S_load=S_load+P_loadjj+1i*Q_loadjj;
                S_gen=S_gen+P_genjj+1i*Q_genjj;
            end

            I_ii_all_mag=abs(conj(S_ii_all/obj.results.V(node_ii)));
            Z_agg=Losses_all/(I_ii_all_mag^2);

            x=10;

        end

        function obj=Gen_aggregation_fin(obj)
            Agg_nodes=[3];
            node_ii=2; % 
            P_loadii=obj.Tables.T_LOAD.P(obj.Tables.T_LOAD.bus==node_ii);
            Q_loadii=obj.Tables.T_LOAD.Q(obj.Tables.T_LOAD.bus==node_ii);
            S_load=P_loadii+1i*Q_loadii;
            S_gen=0;
            S_ii_all=0;
            Losses_all=0;
            for jj=1:length(Agg_nodes)
                node_jj=Agg_nodes(jj);
                P_loadjj=obj.Tables.T_LOAD.P(obj.Tables.T_LOAD.bus==node_jj);
                Q_loadjj=obj.Tables.T_LOAD.Q(obj.Tables.T_LOAD.bus==node_jj);
                P_genjj=obj.Tables.T_VSC.P(obj.Tables.T_VSC.bus==node_jj);
                Q_genjj=obj.Tables.T_VSC.Q(obj.Tables.T_VSC.bus==node_jj);
                S_ii_jj=obj.results.S_tt(node_ii,node_jj);
                S_jj_ii=obj.results.S_tt(node_jj,node_ii);
                Losses_all=Losses_all+S_ii_jj+S_jj_ii;
                S_ii_all=S_ii_all+S_ii_jj;
                S_load=S_load+P_loadjj+1i*Q_loadjj;
                S_gen=S_gen+P_genjj+1i*Q_genjj;
            end


            Load_agg=S_load+Losses_all;

            x=10;

        end
    end
end

