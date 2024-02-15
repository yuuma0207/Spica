%Spica
%���P�b�g�̃p�����[�^�ݒ�y�ьv�Z�p�N���X
%-------------------------------------------------------------------------%
classdef MainSolver
    properties
        %������
        debug
        
        %-------------------Parameters of Rocket---------------------------
        %-----parameter file-----
        param_fn = '';          %�����\�t�@�C����
        param_path = '';        %�����\�C���f�b�N�X
        thrust_fn = '';         %���͗����t�@�C����
        thrust_path = '';       %���͗����C���f�b�N�X
        Xcp_fn = '';            %���͒��S�ʒu�f�[�^�t�@�C����
        dir_scr = '';           %Script�t�H���_
        
        %-----structure-----
        L = 0;                  %�@�̑S��
        Xs = 0;                 %�\���d�S�ʒu
        Ms = 0;                 %�\������
        d = 0;                  %�@�̊O�a
        Xl1 = 0;                %�㕔�����`���O�ʒu
        Xl2 = 0;                %���������`���O�ʒu
        Xl = 0;                 %�ŏ㕔�c�������`���O
        lug_error = 0;          %�����`���O�ʑ��Y���ɂ�郉���`���ɑ΂���X��
        
        Ix = 0;                 %�\���������[�����g(x��)
        Iy = 0;                 %�\���������[�����g(y��)
        Iz = 0;                 %�\���������[�����g(z��)
        py_st = 'p';            %pitch/yaw�Ώ̐�
        
        Xcg = 0;                %�S�@�d�S�ʒu
        M = 0;                  %�S�@����
        I = zeros(3)            %�S�@�����e���\��@�d�S
        
        %-----aerodynamics-----
        Xcp_a = 0;              %���͒��S(�Ό}���p)
        Xcp_b = 0;              %���͒��S(�Ή�����p)
        Xac_a = 0;              %��͒��S�ʒu(�Ό}���p)
        Xac_b = 0;              %��͒��S�ʒu(�Ή�����p)
        P_Fa = 'p';             %��͍�p�_
        Cx_a = 0;               %�ڐ����͌W��(�Ό}���p)
        Cx_b = 0;               %�ڐ����͌W��(�Ή�����p)
        Cz_a = 0;               %�@�����͔��W��
        Cy_b = 0;               %���������͔��W��
        Cac_p = 0;              %�s�b�`���[�����g�W��@��͒��S
        Cac_y = 0;              %���[���[�����g�W��@��͒��S
        Cmq = 0;                %�s�b�`�������[�����g�W��
        Cnr = 0;                %���[�������[�����g�W��
        Clp = 0;                %���[���������[�����g�W��
        
        %-----recovery-----
        pn = 1;                 %�p���V���[�g�̒i��
        n_para = 1;             %n�i��
        Vz_para = 0;            %�p���~�����x
        delay = 0;              %���_����̊J�P�x�ꎞ��
        h_para = 0;             %�J�P���x
        Cd_para = 0;            %�p���R�͌W��
        S_para = 0;             %�p�����e�ʐ�
        D_para = zeros(3,1);    %�p���R��
        l_cord = 1;             %�V���b�N�R�[�h��
        
        %-----engine-----
        Xta = 0;                %�^���N�d�S�ʒu
        Xox = 0;                %�_���܏d�S�ʒu
        Xf = 0;                 %�R���d�S�ʒu
        Xn = 0;                 %�G���W���S��
        Lf = 0;                 %�O���C������
        Ln = 0.03;              %�m�Y������
        df = 0;                 %�R���O�a
        Mf0 = 0;                %�R������
        Mf1 = 0;                %�R�Č�R������
        rho_f = 0;              %�R�����x
        dox = 0;                %�_���܉~�����a
        Lox = 0;                %�_���܍���
        Mta = 0;                %�^���N����
        Mox0 = 0;               %�_���܎���
        rho_ox = 0;             %�_���ܖ��x
        tb = 0;                 %�G���W���R�Ď���
        Mox_d = 0;              %�_���܎��ʗ���
        Mf_d = 0;               %�R�����ʗ���
        Xox0 = 0;
        
        thrust = [];            %���̓f�[�^
        Ts = 0;                 %���̓f�[�^�T�C�Y
        T = zeros(3,1);         %���̓x�N�g��
        T_set = 'Yes';          %���͗����̓Ǎ��̗L��
        T_euler = zeros(3,1);   %�G���W���̃~�X�A���C�������g(�@�̍��W�n�ɑ΂���I�C���[�p�\��)
        
        %-----multi-stage-----
        N = 1;                  %�S�i��
        Ns = 1;                 %�X�e�[�W�ԍ�
        type = 'r';             %�X�e�[�W���
        mode_burn = 't';        %�_�Ώ���
        h_burn = 0;             %�_�΍��x
        t_burn = 0;             %�_�Ύ���
        mode_sep = 't';         %���X�e�[�W��������
        Mn = 0;                 %���X�e�[�W����
        h_sep = 0;              %��i�������x
        t_sep = 0;              %��i��������
        
        %-----wind model-----
        tmp = 0;                %�C��
        Pa = 0;                 %��C��
        rho_a = 0;              %��C���x
        wind_model = 'PowerLaw';   %�������f���I��
        Z0 = 1;                 %�����������荂�x
        n = 1;                  %�������z�W��
        tl = 0;                 %�ŏ㎞��
        
        %-----MSM(NetCDF�`��)�֘A-----
        MSM_fn = '';            %MSM�f�[�^�t�@�C����
        MSM_lon = [];           %�o�x
        MSM_lat = [];           %�ܓx
        MSM_p = [];             %�C��
        MSM_time = [];          %����
        MSM_z = [];             %�W�I�|�e���V�������x
    	MSM_u = [];             %xE��������
        MSM_v = [];             %yE��������
        MSM_temp = [];          %�C��
        
        %-----launch site-----
        L_l = 0;                %�����`����
        z_l = 0;                %�����`����[���x
        L_bottom = 0;           %���˔Ƌ@�̌�[�̋���
        h_plate = 1;            %���˔̒n�ォ��̍���
        lat = 0;                %�˓_�ܓx
        lon = 0;                %�˓_�o�x
        g1 = 9.8;               %�ˏ�̒n�\�ʂł̏d�͉����x
        mgd = 0;                %���C�Ίp
        mode_azm = 'Magnetic';  %�����
        Amgd = zeros(3);        %�^���������ϊ��s��(���ΐ�)
        land_h = 0;             %���n���x
        
        %-----condtion-----
        mode_landing = 'Hard';      %�~�����[�h
        descent_model = 'Vw_model'; %�����������f��
        azm = 0;                    %�ŏ���ʊp
        elev = 0;                   %�ˊp
        q_l = zeros(4,1);           %�����`���p��
        Vw0 = 0;                    %�����
        Wpsi = 0;                   %����
        Vw_s = zeros(3,1);          %������x�N�g��
        
        %-----constants-----
        g0 = 9.80655;           %�d�͉����x
        gasR = 287.05287;       %�C�̒萔[Nm/KgK]
        gamma = 1.403;          %��C�̔�M��
        
        %-----geographi coordinates-----
        Re = 6378137;           %�n���ԓ����a
        f = 1/298.257222101;    %�n���̝G����
        e = 0;                  %�ȉ~�f�ʂ̗��S��
        lat_d = 0;              %�ܓx1�x�ɑΉ����鋗��
        lon_d = 0;              %�o�x1�x�ɑΉ����鋗��
        
        %-----free variables-----
        %�K�v�ȕϐ����g�p�҂����̓s�x�t�������Ă�������
        theta_st = 0;           %�p���p�̏�ԕϐ�
        t_theta = 0;            %���������Ƃ̎p���p��30deg�𒴂��鎞��
        theta_top = 0;          %���_�̎p���p
        z_sep = 0;              %��i�������x
        Cx2 = 0;                %�ڐ����͌W��(������)
        kz2 = 0;                %�@�����͔��W��(������)
        Xcp2 = 0;               %���͒��S�ʒu(������)
        alpha_max = 0;          %�ő�}���p
        beta_max = 0;           %�ő剡����p
        z_ns = 0;               %��i�������x
        theta_ns = 0;           %��i�������p���p
        alpha_lc = 0;           %�����`�N���A���}���p
        beta_lc = 0;            %�����`�N���A��������p
        AngleAccel_max = 0;     %�ő�p�����x(pitch,yaw�����̍��v)
        
        %--------------------------Variables-------------------------------
        %-----times-----
        t = 0;                  %����
        t0 = 0;                 %��������
        t_max = 100;            %�v�Z�ł��؂莞��
        freq = 200;             %���͗����̃T���v�����O���[�g
        delta_t = 0.005;        %���ԍ��ݕ�
        
        %-----flight status-----
        Mox_st = 0;             %�_���ܗ��o���
        Mf_st = 0;              %�R�����o���
        burn_st = 0;            %�R�ď��
        launch_clear_st = 0;    %�����`�N���A���
        top_st = 0;             %���_���B���
        para_st = 0;            %�p���J�����
        sep_st = 0;             %���i�������
        accel_st = 0;           %�����x�̏����̕s�A�����ɑ΂����ԕϐ�
        
        %-----variables-----
        xe0 = zeros(3,1);       %�����ʒu@�n��n
        Ve0 = zeros(3,1);       %�������x@�n��n
        q0 = zeros(4,1);        %�����p��(�N�H�[�^�j�I��)
        omega0 = zeros(3,1);    %�����p���x
        Mf = 0;                 %�R������
        Mox = 0;                %�_���܎���
        
        xe = zeros(3,1);        %�ʒu�x�N�g��@�n��n
        Ve = zeros(3,1);        %��Α��x�x�N�g��@�n��n
        q = zeros(4,1);         %�N�H�[�^�j�I��
        omega = zeros(3,1);     %�p���x�x�N�g��
        
        xe_dot = zeros(3,1);
        Ve_dot = zeros(3,1);
        q_dot = zeros(4,1);
        omega_dot = zeros(3,1);
        M_dot = 0;
        I_dot = zeros(3);
        
        Vw = zeros(3,1);        %�����x�N�g��
        Va = zeros(3,1);        %�΋C���x�x�N�g��
        Van = 0;                %�΋C���x�̑傫��
        Vs = 0;                 %����
        Mach = 0;               %�}�b�n��
        dp = 0;
        
        %-----aerodynamic forces-----
        Xcp = zeros(3,1);       %���͒��S�ʒu
        a = 0;                  %���͒��S�ʒu�f�[�^(�Ό}���p)
        b = 0;                  %���͒��S�ʒu�f�[�^(�Ή�����p)
        Xac = zeros(3,1);       %��͒��S�ʒu
        Cf = zeros(3,1);        %��͌W��
        Cac = zeros(3,1);       %��̓��[�����g�W��@��͒��S
        Cmd = zeros(3,1);       %��͌������[�����g�W��
        Fa = 0;                 %��C�̓x�N�g��
        Ma = 0;                 %��̓��[�����g�x�N�g��
        Mj = 0;                 %�W�F�b�g�_���s���O���[�����g�x�N�g��
        
        %-----Feature Values-----
        Vn_lc = 0;              %�����`�N���A���x
        t_lc = zeros(1,2);      %�����`�N���A����
        Van_max = 0;            %�ő�΋C���x
        Mach_max = 0;           %�ő�}�b�n��
        Accel_max = 0;          %�ő�����x
        top = zeros(3,1);       %���_���W
        z_max = 0;              %�ō����x
        t_top = 0;              %���_���B����
        Van_top = 0;            %���_�΋C���x
        t_para = 0;             %�J�P����
        Va_para = 0;            %�J�P���΋C���x
        t_landing = 0;          %���n����
        dp_max = 0;
        dp_max_t = 0;
        dp_max_z = 0;
        
        
        para_t = 9;
        
        %-----------------------------Solver-------------------------------
        %-----flags-----
        %�\���o�[����
        ode_flag = [1 1];
        definitive_value = true;
        
        %-----ode_adams�̏o��-----
        ts = [];
        xs = [];
        ys = [];
        fs = [];
        
        log_tmp = [];
        
    end
    
    methods
        function obj = MainSolver(gs, stage)    %�e��p�����[�^�t�@�C���Ǎ�, �����l�E�萔�ݒ�
            %GS�N���X�̃v���p�e�B����
            gs_list = properties(gs);
            ms_list = properties(obj);
            list = ismember(ms_list,gs_list);
            list_true = ms_list(list==1);
            for i = 1:size(list_true,1)
                obj.(list_true{i,:}) = gs.(list_true{i,:});
                if size(obj.(list_true{i,:}),1) > 1
                    obj.(list_true{i,:}) = obj.(list_true{i,:})(stage,:);
                end
            end
            
            %�p�����[�^�Ǎ�
            cd(obj.param_path)
            if isfile(obj.param_fn)
                %���C���̕\��
                param_raw = readcell(obj.param_fn,...
                    'Sheet','Spica','Range','A1:B66');
                for k = 1:size(param_raw, 1)
                    if (ischar(param_raw{k, 1})) && (param_raw{k, 1} ~= "")
                        obj.(param_raw{k, 1}) = param_raw{k, 2};
                    end
                end
                
                if obj.type == 'r'
                    switch obj.py_st
                        case 'p'        %pitch�n�̃p�����[�^���g�p
                            obj.Xcp = [0; obj.Xcp_a; obj.Xcp_a];
                            if strcmp(obj.P_Fa,'d')
                                obj.a = readmatrix(obj.Xcp_fn);
                                a0 = obj.a(:,ismember(obj.a(1,:),0));
                                a1 = obj.a(:,ismember(obj.a(1,:),0)==0);
                                if isempty(a0)
                                    obj.a = [-fliplr(a1(1,:)), a1(1,:);
                                        fliplr(a1(2,:)), a1(2,:)];
                                else
                                    obj.a = [-fliplr(a1(1,:)), a0(1,:), a1(1,:);
                                        fliplr(a1(2,:)), a0(2,:), a1(2,:)];
                                end
                            else
                                obj.a = obj.Xcp;
                            end
                            obj.b = obj.a;
                            obj.Xac = [0; obj.Xac_a; obj.Xac_a];
                            obj.Cf = [obj.Cx_a; obj.Cz_a; obj.Cz_a];
                            obj.Cac = [0; obj.Cac_p; obj.Cac_p];
                            obj.Cmd = [obj.Clp; obj.Cmq; obj.Cmq];
                        case 'y'        %yaw�n
                            obj.Xcp = [0; obj.Xcp_b; obj.Xcp_b];
                            if strcmp(obj.P_Fa,'d')

                                obj.b = readmatrix(obj.Xcp_fn);
                                b0 = obj.b(:,ismember(obj.b(1,:),0));
                                b1 = obj.b(:,ismember(obj.b(1,:),0)==0);
                                if isempty(b0)
                                    obj.b = [-fliplr(b1(1,:)), b1(1,:);
                                        fliplr(b1(2,:)), b1(2,:)];
                                else
                                    obj.b = [-fliplr(b1(1,:)), b0(1,:), b1(1,:);
                                        fliplr(b1(2,:)), b0(2,:), b1(2,:)];
                                end
                            else
                                obj.b = obj.Xcp;
                            end
                            obj.a = obj.b;
                            obj.Xac = [0; obj.Xac_b; obj.Xac_b];
                            obj.Cf = [obj.Cx_a; obj.Cy_b; obj.Cy_b];
                            obj.Cac = [0; obj.Cac_y; obj.Cac_y];
                            obj.Cmd = [obj.Clp; obj.Cnr; obj.Cnr];
                        case 'n'        %�Ώ̐��Ȃ�
                            obj.Xcp = [0; obj.Xcp_a; obj.Xcp_b];
                            if strcmp(obj.P_Fa,'d')
                                obj.a = readmatrix(obj.Xcp_fn);
                                a0 = obj.a(:,ismember(obj.a(1,:),0));
                                a1 = obj.a(:,ismember(obj.a(1,:),0)==0);
                                if isempty(a0)
                                    obj.a = [-fliplr(a1(1,:)), a1(1,:);
                                        fliplr(a1(2,:)), a1(2,:)];
                                else
                                    obj.a = [-fliplr(a1(1,:)), a0(1,:), a1(1,:);
                                        fliplr(a1(2,:)), a0(2,:), a1(2,:)];
                                end
                                obj.a = Xcp;
                            end
                            if strcmp(obj.P_Fa,'d')
                                cd(obj.param_path)
                                [obj.Xcp_fn,path] = uigetfile('*.xlsx; *.csv; *.txt','Xcp_data vs beta');
                                obj.Xcp_fn = char(obj.Xcp_fn);
                                cd(path)
                                obj.b = readmatrix(obj.Xcp_fn);
                                b0 = obj.b(:,ismember(obj.b(1,:),0));
                                b1 = obj.b(:,ismember(obj.b(1,:),0)==0);
                                if isempty(b0)
                                    obj.b = [-fliplr(b1(1,:)), b1(1,:);
                                        fliplr(b1(2,:)), b1(2,:)];
                                else
                                    obj.b = [-fliplr(b1(1,:)), b0(1,:), b1(1,:);
                                        fliplr(b1(2,:)), b0(2,:), b1(2,:)];
                                end
                            else
                                obj.b = Xcp;
                            end
                            obj.Xac = [0; obj.Xac_a; obj.Xac_b];
                            obj.Cf = [obj.Cx_a; obj.Cy_b; obj.Cz_a];
                            obj.Cac = [0; obj.Cac_p; obj.Cac_y];
                            obj.Cmd = [obj.Clp; obj.Cmq; obj.Cnr];
                    end
                end
                
                obj.e = obj.f*(2-obj.f);
                R = obj.Re*(1-obj.f)^2/(1+obj.e*cosd(obj.lat));
                obj.lat_d = R*pi/180/3600;
                obj.lon_d = R*cosd(obj.lat)*pi/180/3600;
                
                %�v�Z�̈��
                obj.I = readmatrix(obj.param_fn,...
                    'Sheet', 'Spica', 'Range', 'E1:G3');    %�����e���\��
                
                %�p���V���[�g
                para_raw = readcell(obj.param_fn,...
                    'Sheet', 'Spica', 'range', 'D5:I10');
                for k = 1:6
                    if (ischar(para_raw{k, 1})) && (para_raw{k, 1} ~= "")
                        obj.(para_raw{k, 1}) = cell2mat(para_raw(k, 2:obj.pn+1));
                    end
                end
            else
                warning('Parameter File is NOT FOUND!')
            end
            
            obj.delta_t = 1./obj.freq;           %�v�Z��
            
            switch obj.type
                case 'r'
                    obj.M = obj.Ms + obj.Mf0 + obj.Mox0 + obj.Mta;  %�@�̎���
                    obj.Xcg = (obj.Xs*obj.Ms + obj.Xf*obj.Mf0 + obj.Xox*obj.Mox0 + obj.Xta*obj.Mta)/obj.M;  %�@�̏d�S
                case 'p'
                    obj.M = obj.Ms;
                    obj.Xcg = obj.Xs;
            end
            
            %���̓f�[�^�Ǎ�
            switch obj.type
                case 'r'
                    if strcmp(obj.T_set, 'Yes')
                        cd(obj.thrust_path)
                        obj.thrust = readmatrix(obj.thrust_fn);
                        obj.Ts = size(obj.thrust,1);         %���̓f�[�^�̃T�C�Y
                        if obj.Ns == 1
                            obj.t_burn = 0;
                            obj.burn_st = 1;
                        end
                    end
                case 'p'
                    obj.thrust = 0;
                    obj.t_burn = 0;
            end
            
            %���i��
            if obj.Ns < obj.N
                obj.sep_st = 1;
            end
            cd(obj.dir_scr)
        end
        
        function obj = simulation(obj,Xin)  %�V�~�����[�V�����̃��C���֐�
            switch obj.mode_azm
                case 'Magnetic'
                    obj.Amgd = [cosd(obj.mgd) sind(obj.mgd) 0;
                        -sind(obj.mgd) cosd(obj.mgd) 0;
                        0 0 1];
                case 'Due'
                    obj.Amgd = eye(3);
            end
            
            switch obj.wind_model
                case 'PowerLaw'
                    obj.Vw_s = -[obj.Vw0*cosd(obj.Wpsi); obj.Vw0*sind(obj.Wpsi); 0];
                case 'MSM'
                    %MSM�f�[�^�̓Ǎ�
                    %���s��w�������������̃f�[�^�x�[�X(�ȉ���URL)����f�[�^��\�ߎ擾
                    %http://database.rish.kyoto-u.ac.jp/arch/jmadata/gpv-netcdf.html
                    obj.MSM_fn = uigetfile('*.nc','MSM_data');
                    obj.MSM_fn = char(obj.MSM_fn);
                    MSM_name = {'lon','lat','p','time','z','u','v','temp'};
                    for k = 1:size(MSM_name,2)
                        obj.(strcat('MSM_',MSM_name{k})) = ncread(obj.MSM_fn, MSM_name{k});
                    end
                otherwise
                    warning('wind_model is NOT SELECTED!')
            end
            
            %������
            %-----�t���O-----
            obj.ode_flag = [1 1];
            %ode_flag(1)��0�̊�ode_adams�����Z(t=t0+t_max�Ń��~�b�g)
            %ode_flag(2)��0��Adams�@�̎���k�����Z�b�g
            
            %-----��ԕϐ�-----
            %�������������s�A���ɂȂ�Ƃ��Ȃǂɏ�Ԃ��C���N�������g����
            %����
            obj.Mox_st = 0;
            obj.Mf_st = 0;
            
            %�p���J�P���
            obj.para_st = 0;
            
            %-----�������-----
            %����
            obj.t0 = Xin(14);
            
            %��ԕϐ�
            obj.launch_clear_st = 0;    %�����`�N���A���
            obj.Xl = obj.Xl1;           %�����`���O�ʒu
            obj.top_st = 0;             %���_���B���
            
            %�ŏ����
            obj.elev = -abs(obj.elev);      %�G���[����̂��ߋ����I�ɕ��ɏ㏑��
            if obj.Ns == 1            %�����`�N���A���x
                obj.z_l = obj.L_l * sind(-obj.elev) + obj.h_plate;
            else
                obj.L_l = obj.L - obj.Xl1 + obj.h_plate;
            end
            
            %�p��
            if obj.Ns == 1
                %�����`���p��
                obj.q_l = [cosd(obj.elev/2)*cosd(obj.azm/2);
                           -sind(obj.elev/2)*sind(obj.azm/2);
                           sind(obj.elev/2)*cosd(obj.azm/2);
                           cosd(obj.elev/2)*sind(obj.azm/2)];
                
                %�����@�̎p��
                obj.q0 = [cosd(obj.lug_error/2)*cosd(obj.elev/2)*cosd(obj.azm/2)+sind(obj.lug_error/2)*sind(obj.elev/2)*sind(obj.azm/2);
                          sind(obj.lug_error/2)*cosd(obj.elev/2)*cosd(obj.azm/2)-cosd(obj.lug_error/2)*sind(obj.elev/2)*sind(obj.azm/2);
                          cosd(obj.lug_error/2)*sind(obj.elev/2)*cosd(obj.azm/2)+sind(obj.lug_error/2)*cosd(obj.elev/2)*sind(obj.azm/2);
                          cosd(obj.lug_error/2)*cosd(obj.elev/2)*sind(obj.azm/2)-sind(obj.lug_error/2)*sind(obj.elev/2)*cosd(obj.azm/2)];
                
            else
                obj.q0 = Xin(7:10);
            end
            
            obj.Xox0 = obj.Xox;

            %�ʒu
            l = obj.L - obj.Xcg + obj.L_bottom;             %���˔Ƌ@�̏d�S�̋���
            obj.xe0 = quaternion.q_rot([l;0;0],obj.q0)+Xin(1:3)+[0;0;obj.h_plate];
            
            %���͕Ό�
            obj.T_euler = [0; 0; 0];
            
            %-----�ܓx�E�o�x�֘A-----
            
            %���x,�p���x
            obj.Ve0 = Xin(4:6);
            obj.omega0 = Xin(11:13);
            
            %�����l
            obj.Vn_lc = 0;          %�����`�N���A���x
            obj.t_lc = zeros(1,2);  %�����`�N���A����
            obj.Van_max = 0;        %�ő�΋C���x
            obj.Mach_max = 0;       %�ő�}�b�n��
            obj.Accel_max = 0;      %�ő�����x
            obj.top = zeros(3,1);   %���_���W
            obj.t_top = 0;          %���_���B����
            obj.Van_top = 0;        %���_�΋C���x
            obj.Va_para = 0;        %�J�P���΋C���x
            
            %���R�ϐ��L�q��
            obj.theta_st = 0;
            obj.z_sep = 0;
            obj.Cx2 = 1.23738;
            obj.kz2 = 8.40665;
            obj.Xcp2 = 1.51571;
            obj.alpha_max = 0;
            obj.beta_max = 0;
            obj.AngleAccel_max = 0;
            
            switch obj.type
                case 'r'
                    X0 = [obj.xe0; obj.Ve0; obj.q0; obj.omega0; obj.Mox0; obj.Mf0];
                case 'p'
                    X0 = [obj.xe0; obj.Ve0; obj.q0; obj.omega0; 0; 0];
            end
            
            [ts, xs, ys, fs, obj] = ode_adams(obj, X0, obj.t0, obj.t_max, obj.delta_t);
            obj.ts = ts;
            obj.xs = xs;
            obj.ys = ys;
            obj.fs = fs;
            
        end
        
        function [dx x y ode_flag, obj] = dynamics(obj, t, X, definitive_value)
            %-----�ϐ�------
            obj.t = t;
            obj.xe = X(1:3);
            obj.Ve = X(4:6);
            obj.q = X(7:10);
            obj.omega = X(11:13);
            obj.Mox = X(14);
            obj.Mf = X(15);
            obj.definitive_value = definitive_value;

            %-----�����-----
            if obj.xe(3)-obj.land_h > 0
                z = obj.xe(3);
            else
               z = 0;
            end

            [obj.tmp, obj.Vs, obj.Pa, obj.rho_a] = atmos(z);
            g = obj.g0 * (obj.Re/(obj.Re + z))^2;         %���xz�ł̏d�͉����x

            switch obj.wind_model
                case 'PowerLaw'
                    obj.Vw = (z/obj.Z0)^(1/obj.n) * obj.Vw_s;
                case 'MSM'
                    lon_num = interp1(obj.MSM_lon, 1:size(obj.MSM_lon,1), obj.lon, 'nearest');
                    %obj.longitude�ɑΉ������x�N�g��NetCDF_lon�̃C���f�b�N�X�𓾂�
                    lat_num = interp1(obj.MSM_lat, 1:size(obj.MSM_lat,1), obj.lat, 'nearest');
                    %obj.latitude�ɑΉ������x�N�g��NetCDF_lat�̃C���f�b�N�X�𓾂�
                    time_num = interp1(obj.MSM_time, 1:size(obj.MSM_time,1), obj.tl, 'nearest');
                    %obj.launch_time�ɑΉ������x�N�g��NetCDF_time�̃C���f�b�N�X�𓾂�

                    MSM_z = obj.MSM_z(lon_num, lat_num, :, time_num);
                    MSM_z = MSM_z(:);
                    MSM_u = obj.MSM_u(lon_num, lat_num, :, time_num);
                    MSM_u = MSM_u(:);
                    MSM_v = obj.MSM_v(lon_num, lat_num, :, time_num);
                    MSM_v = MSM_v(:);
                    MSM_temp = obj.MSM_temp(lon_num, lat_num, :, time_num);
                    MSM_temp = MSM_temp(:);

                    p_num = interp1(MSM_z, 1:size(MSM_z,1), z, 'spline');     %���xZ�ƃW�I�|�e���V�������xNetCDF_z�͕ʂł́H
                        %z�ɑΉ������x�N�g��MSM_p�̃C���f�b�N�X�𓾂�
                        %z�̈�����(lon, lat, p, time)�ɑΉ�����ԍ�������
                    Pa = interp1(1:size(obj.MSM_p,1), obj.MSM_p, p_num, 'spline') * 100; %[hPa]��[Pa]�ɕϊ�
                    u = interp1(1:size(MSM_u,1), MSM_u, p_num, 'spline');
                    v = interp1(1:size(MSM_v,1), MSM_v, p_num, 'spline');
                    obj.tmp = interp1(1:size(MSM_temp,1), MSM_temp, p_num, 'spline');

                    obj.Vw= obj.Amgd * [u; v; 0];

                    obj.Vs = sqrt( obj.tmp * obj.gamma * obj.gasR );
                    obj.rho_a = Pa / obj.gasR / obj.tmp;
            end
           
            obj.ode_flag = [1 0];        %�\���o�[����t���O�̏�����
            
            %-----�p��-----
            obj.q = obj.q/norm(obj.q);
            q_inv = quaternion.q_inv(obj.q);
            body = quaternion.q_rot([1;0;0], obj.q);
            theta = acosd(dot(body,[0;0;1]));
            if obj.theta_st == 0
                if abs(theta) < 30
                    obj.t_theta = obj.t;
                else
                    obj.theta_st = 1;
                end
            end
           
            %-----����-----
            %�_���܎��ʕω�
            if obj.Mox_st == 0
                if obj.Mox - obj.Mox_d*obj.delta_t <= 0
                    obj.ode_flag(2) = 1;
                    obj.Mox_st = 1;
                end
            else
                obj.Mox_d = 0;
                obj.Mox = 0;
            end
            
            %�R�����ʕω�
            if obj.Mf_st == 0
               if obj.Mf - obj.Mf_d*obj.delta_t <= obj.Mf1
                   obj.ode_flag(2) = 1;
                   obj.Mf_st = 1;
               end
            else
               obj.Mf_d = 0;
               obj.Mf = obj.Mf1;
            end
            
            %���i��-����
            if obj.sep_st == 1
                switch obj.mode_sep
                    case 't'
                        if obj.t >= obj.t_sep
                            obj.Ms = obj.Ms - obj.Mn;
                            obj.z_ns = z;
                            obj.theta_ns = theta;
                            %obj.Cf = [obj.Cx2; obj.kz2; obj.kz2];
                            %obj.Xcp = [0; obj.Xcp2; obj.Xcp2];
                            obj.ode_flag(2) = 1;
                        end
                    case 'a'
                        if obj.xe(3) >= obj.h_sep
                            obj.Ms = obj.Ms - obj.Mn;
                            obj.ode_flag(2) = 1;
                        end
                end
            end
            
            obj.M_dot = -(obj.Mox_d + obj.Mf_d);
            obj.M = obj.Ms + obj.Mta + obj.Mox + obj.Mf;
            
            
            %-----�d�S-----
            if obj.dox ~= 0
%                 lox_d = -obj.Mox_d/(obj.rho_ox*pi*obj.dox^2/4);     %�_���܍���
%                 Xox_d = lox_d/2;                                    %�_���܏d�S
                if t < 2.88
                    Xox_d = 0.380/2/2.88;
                else
                    Xox_d = 0;
                end
                obj.Xox = obj.Xox0 + Xox_d*t;              %�_���܏d�S�ʒu
            else
                obj.Xox = 0;
                Xox_d = 0;
            end
            
            %�S�@�d�S�ʒu
            obj.Xcg = (obj.Xox*obj.Mox + obj.Xta*obj.Mta + obj.Xf*obj.Mf + obj.Xs*obj.Ms)/obj.M;
            %�d�S�ʒu�ω�
            Xcg_dot = ((obj.Mf_d*obj.Xf+obj.Mox_d*obj.Xox+obj.Mox*Xox_d))/obj.M;
            
            %-----�����e���\��-----
            %�R���������[�����g
            If = obj.Mf/24 * (2*obj.Lf^2 + 3*(1-2*obj.Mf/(obj.rho_f*pi*obj.Lf^2)) + 24*(obj.Xf-obj.Xcg)^2);
            If_dot = obj.Mf*(obj.Mf_d/(2*pi*obj.rho_f*obj.Lf)+2*(obj.Xf-obj.Xcg)*(-Xcg_dot))-obj.Mf_d*((2*obj.Lf^2+3)/24+(obj.Xf-obj.Xcg)^2);
            
            %�_���܊������[�����g
            if obj.dox ~=0
                Iox = obj.Mox/48 * (64+(obj.Mox/obj.rho_ox*pi*obj.dox^2)^2 + 3*obj.dox^2 + 48*(obj.Xox-obj.Xcg)^2);
                Iox_dot = 2*obj.Mox*((obj.Xox-obj.Xcg)*(Xox_d-Xcg_dot)-2*obj.Mox_d/(obj.rho_ox*pi*obj.dox^2)^2)-obj.Mox_d*(obj.dox^2/16+(obj.Xox-obj.Xcg)^2);
            else
                Iox = 0;
                Iox_dot = 0;
            end
                
            %�S�@�����e���\��
            obj.I = [obj.Ix, 0, 0;
                0, obj.Iy+If+Iox, 0;
                0, 0, obj.Iz+If+Iox];
            
            %�����e���\���ω�
            obj.I_dot = [0, 0, 0;
                0, If_dot+Iox_dot, 0;
                0, 0, If_dot+Iox_dot];
            
            %-----����-----
            switch obj.type
                case 'r'
                    if obj.burn_st == 0
                        if strcmp(obj.mode_burn, 'h') && (z > obj.h_burn)
                            obj.t_burn = obj.t;
                        else
                            obj.burn_st = 1;
                        end
                    end
                        
                    thrust_num = round((obj.t - obj.t_burn)*obj.freq + 1);
                    if (thrust_num > 0) && (thrust_num <= obj.Ts)
                        Tb = [obj.thrust(thrust_num); 0; 0];
                        Tb = quaternion.q_rot(Tb, quaternion.euler_q(obj.T_euler));
                    else
                        Tb = zeros(3,1);
                    end
                case 'p'
                    Tb = zeros(3,1);
            end

            %����@�n��n
            obj.T = quaternion.q_rot(Tb, obj.q);

            %-----�d��-----
            Ge = [0; 0; -g];                    %�d��@�n��n
            Gb = quaternion.q_rot(Ge, q_inv);   %�d��@�@�̌n
            
            %-----��C��-----
            Vae = obj.Ve - obj.Vw;                                  %�΋C���x@�n��n
            obj.Va = quaternion.q_rot(Vae, q_inv);                  %�΋C���x@�@�̌n
            obj.Van = norm(obj.Va);
            obj.Mach = obj.Van/obj.Vs;                              %�}�b�n��
            S = pi * obj.d^2/4;                                     %��\�ʐ�
            if obj.Van > 0
                alpha = atan2(obj.Va(3), obj.Va(1));                %�}���p
                beta = asin(obj.Va(2)/obj.Van);                     %������p
            else
                alpha = 0;
                beta = 0;
            end
            Fab = - 1/2 * obj.rho_a * S * obj.Van^2 * obj.Cf.*[1; beta; alpha]; %��C����@�@�̌n
            obj.Fa = quaternion.q_rot(Fab, obj.q);                  %��C����@�n��n
            
            %-----��͍�p�_-----
            switch obj.P_Fa
                case 'c'        %���͒��S(���)�g�p
                    r = obj.Xcp;
                    obj.Cac = zeros(3,1);
                case 'd'        %���͒��S(��)�g�p
                    if obj.launch_clear_st == 0
                        r = [0;
                            interp1(obj.a(1,:), obj.a(2,:), 0,'linear','extrap');
                            interp1(obj.b(1,:), obj.b(2,:), 0,'linear','extrap')];
                    else
                        r = [0;
                            interp1(obj.a(1,:), obj.a(2,:), alpha*180/pi,'linear','extrap');
                            interp1(obj.b(1,:), obj.b(2,:), beta*180/pi,'linear','extrap')];
                    end
                    obj.Cac = zeros(3,1);
                case 'a'        %��͒��S�g�p
                    r = obj.Xac;
            end
            
            %-----��s���-----
            if obj.definitive_value
                %���n����
                if (obj.launch_clear_st >= 2) && (obj.top_st > 0) && (obj.xe(3) <= obj.land_h)
                    obj.ode_flag(1) = 0;
                    obj.t_landing = obj.t;
                end
                
                %�ő�΋C���x
                if obj.Van_max < obj.Van
                    obj.Van_max = obj.Van;      %�ő�΋C���x
                    obj.Mach_max = obj.Mach;    %�ő�}�b�n��
                end
                
                %���_���B����
                if obj.top_st == 0
                    if obj.z_max <= obj.xe(3)
                        obj.top = obj.xe;           %���_���W
                        obj.z_max = obj.xe(3);      %�ō����x
                        obj.t_top = obj.t;          %���_���B����
                        obj.Van_top = obj.Van;      %���_�΋C���x
                        obj.theta_top = theta;      %���_�p���p
                    else
                        obj.top_st = 1;
                    end
                    
                    %�ő�����x
                    Aen = norm(obj.Ve_dot);
                    if Aen > obj.Accel_max
                        obj.Accel_max = Aen;
                    end
                    
                    %�ő�p�����x
                    omega_dot_n = norm(obj.omega_dot(2:3));
                    if omega_dot_n > obj.AngleAccel_max
                        obj.AngleAccel_max = omega_dot_n;
                    end
                end
            end
            
            %-----�^��������-----
            %���i
            obj.Ve_dot = Ge + (obj.T + obj.Fa)/obj.M;
            Vb_dot = quaternion.q_rot(obj.Ve_dot, q_inv);
            
            xlb = [obj.Xl-obj.Xcg; 0; 0];
            xl = quaternion.q_rot(xlb, obj.q);
            if obj.launch_clear_st <= 1         %�����`�N���A���
                if obj.launch_clear_st == 0
                    if obj.xe(3)-xl < obj.z_l
                        ql_inv = quaternion.q_inv(obj.q_l);
                        Vl_dot = quaternion.q_rot(obj.Ve_dot, ql_inv);
                        Vl_dot = [Vl_dot(1);0;0];
                        obj.Ve_dot = quaternion.q_rot(Vl_dot, obj.q_l);
                        if obj.accel_st == 0
                            if obj.Ns==1 && obj.Ve_dot(3)<0
                                obj.Ve_dot = zeros(3,1);
                            else
                                obj.ode_flag(2) = 1;
                                obj.accel_st = 1;
                            end
                        end
                    else
                        obj.launch_clear_st = 1;
                        obj.Xl = obj.Xl2;
                        obj.t_lc(1,1) = obj.t;
                    end
                else
                    if obj.xe(3)-xl >= obj.z_l
                        if obj.definitive_value
                            obj.ode_flag(2) = 1;
                            obj.launch_clear_st = 2;
                            obj.Vn_lc = norm(obj.Ve);
                            obj.t_lc(1,2) = obj.t;
                            obj.alpha_lc = alpha;
                            obj.beta_lc = beta;
                        end
                    end
                end
            end
            
            %��]
            l = [0;obj.Xcg;obj.Xcg] - r;
            Mfa = l.*Fab;
            Mfa = [Mfa(1);-Mfa(3);Mfa(2)];
            Mac = 1/2 * obj.rho_a * S * obj.Van^2 * obj.Cac * obj.L;
            Md = 1/4*obj.rho_a*S*obj.Van*obj.Cmd.*[obj.d^2;obj.L^2;obj.L^2].*obj.omega;
                %����͌������[�����g�x�N�g��
            obj.Ma = Mfa + Mac + Md;
                %����̓��[�����g�x�N�g��
            re = obj.Xn - obj.Xcg;
            rt = re - obj.Ln;
            obj.Mj = obj.M_dot * (re*rt) .* obj.omega;          %-diag(obj.I)/obj.M
            
                %���W�F�b�g�_���s���O���[�����g�x�N�g��
            obj.Mj(1) = 0;                      %���[���͖���
            gyro = -cross(obj.omega, obj.I*obj.omega);
            
            if obj.launch_clear_st == 0
                obj.Ma = zeros(3,1);
                obj.Mj = zeros(3,1);
                gyro = zeros(3,1);
                obj.omega = zeros(3,1);
                Mg = zeros(3,1);
            else
                Mg = zeros(3,1);
            end
            
            obj.omega_dot = obj.I\(obj.Ma + obj.Mj + gyro - obj.I_dot * obj.omega + Mg);
            
            %-----�p���V���[�g-----
            if obj.top_st == 1
               if strcmp(obj.mode_landing,'Descent')
                   if obj.para_st+1 <= obj.pn
                       if (obj.delay(obj.para_st+1) <= obj.t-obj.t_top) || (obj.h_para(obj.para_st+1) >= obj.xe(3))
                            obj.ode_flag(2) = 1;
                            obj.para_st = obj.para_st + 1;
                            obj.t_para = obj.t;
                        end
                   else
                       obj.para_st = obj.pn;
                   end

                   if obj.para_st > 0
                       switch obj.descent_model
                           case 'Vw_model'
                               obj.Ve = obj.Vw + [0; 0; -obj.Vz_para(obj.para_st)];
                           case 'Dynamics'
                               obj.D_para = - 1/2 * obj.rho_a * obj.S_para * obj.Van^2 * obj.Cd_para * Vae/obj.Van;
                               obj.Ve_dot = obj.D_para/obj.M + Ge;
                               obj.omega = zeros(3,1);
                       end
                   end
               end
            end
            
            %-----�ʒu�E�p��-----
            obj.xe_dot = obj.Ve;
            P = obj.omega(1);
            Q = obj.omega(2);
            R = obj.omega(3);
            obj.q_dot = 1/2 * [0, -P, -Q, -R;
                               P, 0, R, -Q;
                               Q, -R, 0, P;
                               R, Q, -P, 0] * obj.q;
            
            %-----�ܓx�E�o�x�ϊ�-----
            
            
            %-----���i��-----
            if obj.sep_st == 1
                switch obj.mode_sep
                    case 't'
                        if obj.t >= obj.t_sep
                            obj.sep_st = 2;
                            obj.z_sep = obj.xe(3);
                        end
                    case 'a'
                        if obj.xe(3) >= obj.h_sep
                            obj.sep_st = 2;
                        end
                end
            end
            
            %-----���R�ϐ�-----
            if obj.para_st == 0
                obj.Va_para = obj.Van;
            end
            if obj.launch_clear_st == 1
                if abs(alpha) >= obj.alpha_max
                    obj.alpha_max = abs(alpha);
                end
                if abs(beta) >= obj.beta_max
                    obj.beta_max = abs(beta);
                end
            end

            Fst = (r(2)-obj.Xcg(1))/obj.L*100;
            dp = 1/2*obj.rho_a*obj.Van^2/10^3;
            extra = [Fst; dp];
           if definitive_value
               if dp > obj.dp_max && obj.top_st == 0
                   obj.dp_max = dp;
                   obj.dp_max_t = t;
                   obj.dp_max_z = obj.xe(3);
               end
           end
            
            
            %-----�o��-----
            dx = [ obj.xe_dot;
				obj.Ve_dot;
				obj.q_dot;
				obj.omega_dot;
				-obj.Mox_d;
				-obj.Mf_d];
            x = [obj.xe;
                obj.Ve;
                obj.q;
                obj.omega;
                obj.Mox;
                obj.Mf];
            y = [Vb_dot;
                obj.Va;
                obj.Van;
                alpha;
                beta;
                theta;
                obj.Vw;
                obj.D_para;
                Mg;
                r;
                extra];
            
            ode_flag = obj.ode_flag;
            
        end
    
        function [t, x, y, f, obj] = ode_adams(obj, x0, t0, t_max, delta_t)
            %�\���q�C���q�@�iAdams-Bashforth-Moulton PE(CE)2 �@�j�̃\���o�[
            %
            %definitive_value�ɂ���, 
            %��������t�ɂ�����, ���̃\���o�[��ode_function()�𕡐���(3��)�v�Z����. 
            %���̂���, ���O�ɌĂяo���ꂽ�l�Ƃ̔�r�ɂ���ċɑ�l�E�ɏ��l�����߂悤�Ƃ����Ƃ�, �قړ����l���r����, 
            %���̔����Ȓl�̑召�ɂ��, �ɒl�����܂����߂��Ȃ��\��������. 
            %definitive_value�͓�����t�ōŌ��ode_function()���v�Z���鎞��true,
            %���̑��̎���false�����^������悤�ɂȂ��Ă���. ���̏����Q�l�ɏ������򂵂ċɒl�����߂�Ƃ悢. 

            %Adams�@(Bashforth, Moulton)�̌W��
            k_max = 10;
            bashforth_coefficient = {
                flip(1.0) %�I�C���[�@�̌W��
                flip([3.0; -1.0] / 2.0) %k=2
                flip([23.0; -16.0; 5.0] / 12.0)
                flip([55.0; -59.0; 37.0; -9.0] / 24.0)
                flip([1901.0; -2774.0; 2616.0; -1274.0; 251.0] / 720.0)
                flip([4277.0; -7923.0; 9982.0; -7298.0; 2877.0; -475.0] / 1440.0)
                flip([198721.0; -447288.0; 705549.0; -688256.0; 407139.0; -134472.0; 19087.0] / 60480.0)
                flip([434241.0; -1152169.0; 2183877.0; -2664477.0; 2102243.0; -1041723.0; 295767.0; -36799.0] / 120960.0)
                flip([14097247.0; -43125206.0; 95476786.0; -139855262.0; 137968480.0; -91172642.0; 38833486.0; -9664106.0; 1070017.0] / 3628800.0)
                flip([30277247.0; -104995189.0; 265932680.0; -454661776.0; 538363838.0; -444772162.0; 252618224.0; -94307320.0; 20884811.0; -2082753.0] / 7257600.0)
            };
            moulton_coefficient = {
                flip([1.0; 1.0] / 2.0) %�C���I�C���[�@�̌W��
                flip([1.0; 1.0] / 2.0) %k=2
                flip([5.0; 8.0; -1.0] / 12.0)
                flip([9.0; 19.0; -5.0; 1.0] / 24.0)
                flip([251.0; 646.0; -264.0; 106.0; -19.0] / 720.0)
                flip([475.0; 1427.0; -798.0; 482.0; -173.0; 27.0] / 1440.0)
                flip([19087.0; 65112.0; -46461.0; 37504.0; -20211.0; 6312.0; -863.0] / 60480.0)
                flip([36799.0; 139849.0; -121797.0; 123133.0; -88547.0; 41499.0; -11351.0; 1375.0] / 120960.0)
                flip([1070017.0; 4467094.0; -4604594.0; 5595358.0; -5033120.0; 3146338.0; -1291214.0; 312874.0; -33953.0] / 3628800.0)
                flip([2082753.0; 9449717.0; -11271304.0; 16002320.0; -17283646.0; 13510082.0; -7394032.0; 2687864.0; -583435.0; 57281.0] / 7257600.0)
            };

            ns_max = t_max / delta_t; %�X�e�b�v��ns�̍ő�l(�������[�v������邽��)

            %������
            ns = 1; %�X�e�b�v��ns

            x = zeros(size(x0, 1), ns_max); %�ϐ�x
            x(:,1) = x0; %x�̏����l�Ƃ��ė^����ꂽx0����

            [f0, ~, y0, ~, obj] = obj.dynamics(t0, x(:,1), true); %f, y�̏�����, �u~�v�͊֐��o�̖͂��������Ă���
            f = zeros(size(x0, 1), ns_max); %�ϐ�x�̎��Ԕ���f
            f(:,1) = f0;

            y = zeros(size(y0, 1), ns_max); %�ϐ�y
            y(:,1) = y0;

            ode_flag = [1, 1]; %ode_flag�̏�����(���[�v�̎��s, Step��k�̃��Z�b�g)

            while ode_flag(1) && (ns < ns_max) %ode_flag(1)��0����ns<=ns_max�̂Ƃ����[�v�����s����
                %�X�e�b�v��ns, ����t�̍X�V
                ns = ns + 1;
                t = t0 + (ns-1) * delta_t;
                %Adams�@��Step��k�̌���
                if ode_flag(2) %Step��k�̃��Z�b�g�t���O(ode_flag(2)��0)
                    k = 1;
                else
                    if k < k_max
                        k = k + 1;
                    end
                end

                %�\���q( Adams-Bashforth �@)
                x(:,ns) = x(:,ns-1) + delta_t * my_product( f, bashforth_coefficient{k}, ns-1 );
                [f(:,ns), x(:,ns), ~, ~, obj] = obj.dynamics( t, x(:,ns), false );

                %�C���q1( Adams-Moulton �@)
                x(:,ns) = x(:,ns-1) + delta_t * my_product( f, moulton_coefficient{k}, ns );
                [f(:,ns), x(:,ns), ~, ~, obj] = obj.dynamics( t, x(:,ns), false );

                %�C���q2( Adams-Moulton �@)
                x(:,ns) = x(:,ns-1) + delta_t * my_product( f, moulton_coefficient{k}, ns );
                [f(:,ns), x(:,ns), y(:,ns), ode_flag, obj] = obj.dynamics( t, x(:,ns), true );
            end

            %ode_flag(1)=0�őł��؂�ꂽ�Ƃ��A���̗]�v�ȃf�[�^����������
            t = t0 + ( 0:(ns-1) ) * delta_t; %����t
            x = x(:,1:ns);
            y = y(:,1:ns);
            f = f(:,1:ns);

            %%�⏕�֐�
                function an = my_product(f, coefficient, p)
                    %���ς̌v�Z(f��p����coefficient�����̌��k���ė�x�N�g����؂�o���ē��ς��v�Z����B
                    coefficient_dim = size(coefficient,1);
                    an = f(:, (p - coefficient_dim + 1):p) * coefficient;
                end
        end
    
    end

end

