%Spica
%���ʕ\���p�N���X
%-------------------------------------------------------------------------%
classdef DisplayFigure
    properties
        %-----directory-----
        dir_home = '';              %�z�[���f�B���N�g��(Spica�t�H���_)
        dir_res = '';               %Result�t�H���_
        dir_scr = '';               %Script�t�H���_
        dir_form = '';              %Formats�t�H���_
        dir_param = '';             %ParameterFiles�t�H���_
        dir_thrust = '';            %ThrustData�t�H���_
        dir_ls = '';                %LaunchSite�t�H���_
        
        %-----parameter-----
        param_n = 1;                %�S�i��
        
        %-----calculating condition-----
        base_azm = 'ME';            %����� M:Magnetic T:True, E:East N:North ...
        mode_angle = 'CCW';         %���ʊp�̐����� CW:ClockWise, CCW:CounterClockWise
        view_azm = 'Magnetic';      %�������U�}�̊���ʎ��
        mode_calc = 'Single';       %�v�Z���[�h  Single, Multiple, FallPoint
        mode_landing = 'Hard';      %�~�����[�h  Hard, Descent, Both
        descent_model = 'Vw_model'; %�����~�����f��    Vw_model, Dynamics
        elev_set = zeros(3,1);      %�ˊp�ݒ�
        Vw0_set = zeros(3,1);       %������ݒ�
        Wpsi_set = zeros(3,1);      %�����ݒ�
        mgd = 0;                    %�˓_���C�Ίp(����=��)
        
        %-----calculating-----
        freq = 0;                   %�v�Z���[�g
        elev = [];                  %�ˊp���X�g
        Vw0 = [];                   %�������X�g
        Wpsi = [];                  %�������X�g
        Wpsi_res = [];              %�������X�g(���ʗp)
        elev_n = 0;                 %�ˊp������
        Vw0_n = 0;                  %����������
        Wpsi_n = 0;                 %����������
        
        %-----result-----
        result_name = ["Hard";"Descent"];   %���ʔz��
        Hard = [];                          %�e���v�Z����(MainSolver�N���X)
        Descent = [];                       %�����v�Z����(MainSolver�N���X)
        FP_Hard = [];                       %�c���������_
        FP_Descent = [];                    %�����������_
        ll = [];                            %�o�ܓx�ϊ��p�N���X(lon_lat�N���X)
        
        limit_area = [];
        %-----display-----
        list_fig = "";                      %�o�͂���}�̃��X�g
        fig_size = struct(...
            'path',[100, 50, 800, 700],...
            'point',[100, 50, 800, 700]);   %�}�̃T�C�Y
        ax = struct(...
            'path',struct(...
                'range', [-840, 880, -910, 750],...
                'label', ["Magnetic East [m]";
                        "Magnetic North [m]";
                        "Altitude [m]"],...
                'FontSize', 10,...
                'Color', [0, 0, 0],...
                'FontWeight', 'normal'),...
            'point',struct(...
                'range', [-840, 880, -910, 750],...
                'label', ["True East[m]";"True North[m]"],...
                'FontSize', 10,...
                'Color', [0, 0, 0],...
                'FontWeight', 'normal'));   %���֘A
        lgd = struct(...
            'path', struct(...
                'pos', [0.8, 0.8, 0.1 ,0.1],...
                'FontSize', 10,...
                'TextColor', [1, 0.9999, 0.9999],...
                'FontWeight', 'bold'),...
            'point', struct(...
                'pos', [0.8, 0.8, 0.1 ,0.1],...
                'FontSize', 10,...
                'TextColor', [1, 0.9999, 0.9999],...
                'FontWeight', 'bold'));     %�}��֘A
        
        %flightpath
        event_list = "";                    %��s�C�x���g���X�g
        event = "";                         %��s�C�x���g
        
        %fallpoint
        back_pict = struct(...
            'fn', 'Oshima_201903.jpg',...
            'img', [],...
            'pos', [-830, 870; -907, 750]); %�w�i�}�֘A
        marker = struct(...
            'mode', 'shape',...
            'color', [1, 1, 1; 1, 0, 0],...
            'size', 50,...
            'shape', "o",...
            'mfc', 'flat');                 %�}�[�J�[�֘A
        
        %kml
        kml = struct(...
            'path', struct(...
                'Width',15,...
                'Color', [1 1 1]),...
            'point', struct(...
                'Width',15,...
                'Color', [1 1 1]));         %GoogleEarth�pkml�t�@�C���֘A
        kml_str = struct(...
            'str1', ["<?xml version=""1.0"" encoding=""utf-8""?>";
                "<kml xmlns=""http://www.opengis.net/kml/2.2"">";
                "<Document>"],...
            'str2', "<LineStyle>",...
            'str3', ["</LineStyle>";
                "<ColorStyle>"],...
            'str4', ["</ColorStyle>";
                "<Placemark>";
                "<Snippet maxLines=""0""> </Snippet>";
                "<description> </description>";
                "<MultiGeometry>"],...
            'str5',["<LineString>";
                "<altitudeMode>relativeToGround</altitudeMode>";
                "<coordinates>"],...
            'str_multi',["</coordinates>";
                "</LineString>"],...
            'str6', ["</coordinates>";
                "</LineString>";
                "</MultiGeometry>";
                "</Placemark>";
                "</Document>";
                "</kml>"]);
        
    end
    
    methods
        function obj = DisplayFigure(gs, cc)
            %GeneralSetting�N���X�̃v���p�e�B����
            gs_list = properties(gs);
            df_list = properties(obj);
            list = ismember(df_list,gs_list);
            list_true = df_list(list==1);
            for i = 1:size(list_true,1)
                obj.(list_true{i,:}) = gs.(list_true{i,:});
            end
            
            %Calculation�N���X�̃v���p�e�B����
            cc_list = properties(cc);
            list = ismember(df_list,cc_list);
            list_true = df_list(list==1);
            for i = 1:size(list_true,1)
                obj.(list_true{i,:}) = cc.(list_true{i,:});
            end
            
            obj.ax.point.label = ["True East[m]";"True North[m]"];
            %�e��}�̏o��
            path_tf = ismember(["FlightPath";"KML of FlightPath"], obj.list_fig);
            point_tf = ismember(["FallPoint";"KML of FallPoint"], obj.list_fig);
            if sum(path_tf) > 0
                disp("OutPutting FlightPath Figures...")
                obj.output_path;
            end
            if sum(point_tf) > 0
                disp("OutPutting FallPoint Figures...")
                %�����_�z��̎擾
                for i = 1:size(obj.mode_landing,1)
                    obj.(strcat("FP_",obj.mode_landing{i,:})) = cc.feature.(strcat("FP_",obj.mode_landing{i,:}));
                end
                obj.output_point;
            end
        end
        
        function output_path(obj)       %��s�o�H�}�o��
            ax_set = obj.ax.path;
            lgd_set = obj.lgd.path;
            
            %kml�t�@�C���pcolor�ϊ�
            str_color = "ff";
            for i = 1:3
                tmp = dec2base(obj.kml.path.Color(i)*255,16);
                if size(tmp,2) == 1
                    strcat('0',tmp')
                end
                str_color = strcat(str_color,lower(tmp));
            end
            
            for h = 1:obj.elev_n
                for i = 1:obj.Vw0_n
                    %���ʃt�H���_����
                    if isfolder(obj.dir_res)==0
                        cd(obj.dir_home);
                        mkdir 'Result'
                        obj.dir_res = strcat(obj.dir_home,'/Result');
                    end
                    cd(obj.dir_res)
                    st_elev = strcat(num2str(obj.elev(h)),'deg');
                    if isfolder(st_elev)==0
                        mkdir(st_elev)
                    end
                    cd(st_elev)
                    st_Vw0 = strcat(num2str(obj.Vw0(i)),'m');
                    if isfolder(st_Vw0)==0
                        mkdir(st_Vw0)
                    end
                    cd(st_Vw0)
                    dir_name = cd(obj.dir_scr);
                    
                    for j = 1:obj.Wpsi_n
                        for k = 1:size(obj.mode_landing,1)
                            %�t�@�C����
                            cond_name = strcat("FlightPath_",...
                                num2str(obj.elev(h)),"deg_",...
                                num2str(obj.Vw0(i)),"m_",...
                                num2str(obj.Wpsi(j)),"deg_",...
                                num2str(obj.mode_landing(k)));
                            
                            for l = 1:obj.param_n
                                tmp_res = obj.(obj.mode_landing{k,:})(h,i,j,l);
                                
                                %���W�̎��ԃ��O�𒊏o
                                pos = tmp_res.xs(1:3,:)';
                                
                                %�X�e�[�W�ԍ�
                                str_l = char(num2str(l));
                                switch str_l(1,size(str_l,2))
                                    case '1'
                                        stage_str = strcat(str_l,'st');
                                    case '2'
                                        stage_str = strcat(str_l,'nd');
                                    case '3'
                                        stage_str = strcat(str_l,'rd');
                                    otherwise
                                        stage_str = strcat(str_l,'th');
                                end
                                stage_str = strcat(stage_str, " stage");
                                
                                %-----jpg�t�@�C��-----
                                if ismember("FlightPath", obj.list_fig)
                                    %�v���b�g
                                    plot3(pos(:,1),pos(:,2),pos(:,3),...
                                        'DisplayName',strcat("Ns=",num2str(l)))
                                    hold on
                                    grid on
                                    if l == 1
                                        %�}�̃v���p�e�B���擾
                                        fig_name = strcat('FlightPath_',obj.mode_landing(k,:),...
                                            '(elev=',num2str(obj.elev(h)),',Ns=',num2str(l),...
                                            ',Vw0=',num2str(obj.Vw0(i)),...
                                            ',Wpsi=',num2str(obj.Wpsi(j)),')');
                                        fig = gcf;
                                        fig.Name = fig_name;
                                        fig.Position = obj.fig_size.path;
                                    end
                                    if l == obj.param_n
                                        %���ݒ�
                                        tmp_ax = gca;
                                        xlabel(ax_set.label(1,:),'FontSize',ax_set.FontSize)
                                        ylabel(ax_set.label(2,:),'FontSize',ax_set.FontSize)
                                        zlabel(ax_set.label(3,:),'FontSize',ax_set.FontSize)
                                        tmp_ax.FontSize = ax_set.FontSize;
                                        tmp_ax.XColor = ax_set.Color;
                                        tmp_ax.YColor = ax_set.Color;
                                        tmp_ax.FontWeight = ax_set.FontWeight;
                                        %�}��쐬�E�ݒ�
                                        tmp_lgd = legend;
                                        tmp_lgd.Position = lgd_set.pos;
                                        tmp_lgd.FontSize = lgd_set.FontSize;
                                        tmp_lgd.TextColor = lgd_set.TextColor;
                                        tmp_lgd.FontWeight = lgd_set.FontWeight;
                                    end
                                    eve_list = string(fieldnames(obj.event));
                                    for m = 1:size(eve_list,1)
                                        eve_name = eve_list(m,:);
                                        tmp_eve = obj.event.(eve_name);
                                        if strcmp(tmp_eve.disp,"Yes")
                                            if sum(ismember([obj.mode_landing(k,:),"Both"],tmp_eve.mode_land)) > 0
                                                if tmp_eve.stage==l
                                                    t_eve = tmp_res.(tmp_eve.ref)(1)
                                                    pos_eve = pos(round(t_eve*obj.freq),:);
                                                    text(pos_eve(1),pos_eve(2),pos_eve(3),...
                                                        strcat("\leftarrow",stage_str," ",eve_name))
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                %-----kml�t�@�C��-----
                                if ismember("KML of FlightPath", obj.list_fig)
                                    %�ܓx�o�x�ϊ�
                                    [x_geo, ~] = obj.ll.Vincenty_direct(pos(:,1:2));
                                    x_geo = [x_geo(:,2),x_geo(:,1),pos(:,3)];
                                    if l == 1
                                        kml_fn = strcat(dir_name,"/",cond_name,".kml");
                                        f = fopen(kml_fn, 'w');
                                        write_str = [obj.kml_str.str1;
                                            strcat("<name>",cond_name,"</name>");
                                            obj.kml_str.str2;
                                            num2str(obj.kml.path.Width);
                                            obj.kml_str.str3
                                            str_color;
                                            obj.kml_str.str4];
                                        fprintf(f, '%s\n', write_str);
                                        fclose(f);
                                        f = fopen(kml_fn,'a');
                                    else
                                        fprintf(f, '%s\n', obj.kml_str.str_multi);
                                    end
                                    write_str = [strcat("<name>",stage_str,"</name>");
                                        obj.kml_str.str5];
                                    fprintf(f, '%s\n', write_str);
                                    fprintf(f, '%.12f,%.12f,%.12f\n', x_geo');
                                end
                            end
                            %�}�̕ۑ�
                            if ismember("FlightPath", obj.list_fig)
                                hold off
                                fig_fn = strcat(dir_name,'/',cond_name,'.jpg');
                                saveas(fig, fig_fn)
                                close(fig)
                            end
                            
                            %kml�t�@�C���ۑ�
                            if ismember("KML of FlightPath", obj.list_fig)
                                fprintf(f, '%s\n', obj.kml_str.str6);
                                fclose(f);
                            end
                        end
                    end
                end
            end
        end
        
        function output_point(obj)      %�������U�}�o��
            %�o�͐�̃t�H���_�쐬
            if isfolder(obj.dir_res)==0
                cd(obj.dir_home);
                mkdir 'Result'
                obj.dir_res = strcat(obj.dir_home,"/Result");
            end
            dir_name = repmat(strcat(obj.dir_res,"/Ns_"),obj.param_n,1);
            for i = 1:obj.param_n
                dir_name(i,:) = strcat(dir_name(i,:), num2str(i));
                if isfolder(dir_name(i,:))==0
                    mkdir(strcat(dir_name(i,:)))
                end
            end
            
            %�o�͐ݒ�擾
            ax_set = obj.ax.point;
            lgd_set = obj.lgd.point;
            back_set = obj.back_pict;
            marker_set = obj.marker;
            mkr_list = {'o';'s';'d';'^';'p';'h';'+';'*';'x'}; 
            
            %kml�t�@�C���pcolor�ϊ�
            str_color = "ff";
            for i = 1:3
                tmp = dec2base(obj.kml.point.Color(i),16);
                str_color = strcat(str_color,lower(tmp));
            end
            
            for h = 1:obj.elev_n
                for i = 1:obj.param_n
                    for j = 1:size(obj.mode_landing, 1)
                        str_i = char(num2str(i));
                        switch str_i(1,size(str_i,2))
                            case '1'
                                stage_str = strcat(str_i,'st');
                            case '2'
                                stage_str = strcat(str_i,'nd');
                            case '3'
                                stage_str = strcat(str_i,'rd');
                            otherwise
                                stage_str = strcat(str_i,'th');
                        end
                        cond_name = strcat("FallPoint_",stage_str,"_",...
                            num2str(obj.elev(h)),"deg_",...
                            num2str(obj.mode_landing(j)));
                        
                        if ismember("FallPoint", obj.list_fig)
                            %�}�̍쐬
                            fig_name = strcat('FallPoint_',obj.mode_landing(j,:),...
                                '(elev=',num2str(obj.elev(h)),',Ns=',num2str(i),')');
                            fig = figure('Name',fig_name,'Position',obj.fig_size.point);
                            axis(ax_set.range)
                            hold on
                            image(back_set.pos(1,:),back_set.pos(2,:),back_set.img)
                        end
                        
                        %�v���b�g
                        for k = 1:obj.Vw0_n
                            if ismember("FallPoint", obj.list_fig)
                                %�����_�z��̓Ǎ�
                                fp = obj.(strcat("FP_",obj.mode_landing(j,:)))(:,k,:,i,h);
                                fp = squeeze(fp);
                                if size(fp) == [2,1]
                                    fp = fp';
                                end

                                if strcmp(obj.mode_calc, 'Single')
                                    marker_set.mode =  'shape';
                                end
                                switch marker_set.mode
                                    case 'shape'
                                        mkr = char(mkr_list(k));
                                        scatter(fp(:,1), fp(:,2), marker_set.size,...
                                            marker_set.color(1,:), mkr,...
                                            'LineWidth',1.5,...
                                            'MarkerFaceColor',marker_set.mfc,...
                                            'DisplayName',strcat(num2str(obj.Vw0(k)),'m/s'));
                                    case 'gradation'
                                        color = interp1([1,obj.Vw0_n], marker_set.color,k);
                                        scatter(fp(:,1), fp(:,2), marker_set.size,...
                                            color, marker_set.shape,...
                                            'LineWidth',1.5,...
                                            'MarkerFaceColor',marker_set.mfc,...
                                            'DisplayName',strcat(num2str(obj.Vw0(k)),'m/s'));
                                end

                            end
                            
                            if ismember("KML of FallPoint", obj.list_fig)
                                %-----kml�t�@�C��-----
                                %�ܓx�o�x�ϊ�
                                [x_geo, ~] = obj.ll.Vincenty_direct(fp);
                                x_geo = [x_geo(:,2),x_geo(:,1)];
                                
                                if k == 1
                                    kml_fn = strcat(dir_name(i,:),"/",cond_name,".kml");
                                    f = fopen(kml_fn, 'w');
                                    write_str = [obj.kml_str.str1;
                                        strcat("<name>",cond_name,"</name>");
                                        obj.kml_str.str2;
                                        num2str(obj.kml.point.Width);
                                        obj.kml_str.str3
                                        str_color;
                                        obj.kml_str.str4];
                                    fprintf(f, '%s\n', write_str);
                                    fclose(f);
                                    f = fopen(kml_fn,'a');
                                else
                                    fprintf(f, '%s\n', obj.kml_str.str_multi);
                                end
                                write_str = [strcat("<name>",strcat(num2str(obj.Vw0(k)),"m/s"),"</name>");
                                    obj.kml_str.str5];
                                fprintf(f, '%s\n', write_str);
                                fprintf(f, '%.12f,%.12f\n', x_geo');
                            end
                        end
                        if ismember("FallPoint", obj.list_fig)
                            %���ݒ�
                            tmp_ax = gca;
                            xlabel(ax_set.label(1,:),'FontSize',ax_set.FontSize)
                            ylabel(ax_set.label(2,:),'FontSize',ax_set.FontSize)
                            tmp_ax.FontSize = ax_set.FontSize;
                            tmp_ax.XColor = ax_set.Color;
                            tmp_ax.YColor = ax_set.Color;
                            tmp_ax.FontWeight = ax_set.FontWeight;
                            %�}��쐬�E�ݒ�
                            tmp_lgd = legend('boxoff');
                            tmp_lgd.Position = lgd_set.pos;
                            tmp_lgd.FontSize = lgd_set.FontSize;
                            tmp_lgd.TextColor = lgd_set.TextColor;
                            tmp_lgd.FontWeight = lgd_set.FontWeight;
                            
                            
                            polygon = obj.limit_area.polygon.dist;
                            polygon = [polygon; polygon(1,:)];
                            plot(polygon(:,1),polygon(:,2),'-o','Color','y','MarkerEdgeColor','y',...
                                'MarkerFaceColor','y','DisplayName','Limited Polygon')

                            hold off
                            %�}�̕ۑ�
                            fig_name = strcat(dir_name(i,:),'/',cond_name,'.jpg');
                            saveas(fig, fig_name)
                            close(fig)
                        end
                        
                        if ismember("KML of FallPoint", obj.list_fig)
                            %kml�t�@�C���ۑ�
                            fprintf(f, '%s\n', obj.kml_str.str6);
                            fclose(f);
                        end
                    end
                end
            end
        end
        
        function output_logfig(obj)     %�C�ӕϐ��̃��O��}�ɏo��
            
        end
    end
end