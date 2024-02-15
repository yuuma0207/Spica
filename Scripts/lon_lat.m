classdef lon_lat
   %���ʊp�E�����̌o�ܓx�Ԃ̕ϊ��N���X
   %�n���ȉ~�̂�WGS84(Googl�}�b�v��GPS�Ɠ���)���̗p
   %���ʊp�� �k0deg, CW (����)
   
properties
    a = 6378137.06;         %�n���ԓ����a [m] (WGS84)
    b = 0;                  %�n���ɔ��a [m]
    f = 1/298.257223563;    %�n���̝G���� (WGS84)
    e = 0;                  %�ȉ~�f�ʂ̗��S��
    phi_origin = 0;         %��_�ܓx [deg]
    L_origin = 0;           %��_�o�x [deg]
    U_origin = 0;           %��_�̍X���ܓx [deg]
end

methods
    function obj = lon_lat(pos)     %�R���X�g���N�^���\�b�h
        obj.b = obj.a*(1-obj.f);
        obj.e = sqrt(obj.f*(2-obj.f));
        obj.phi_origin = deg2rad(pos(1));
        obj.L_origin = deg2rad(pos(2));
        obj.U_origin = atan2((1-obj.f) * tan(obj.phi_origin),1);
    end
    
    function [x_geo, alpha_x] = Vincenty_direct(obj, x_dist)
        %��_�̌o�ܓx(obj.phi_origin, obj.L_origin), �ڕW�_���Wx_dist = [x,y]����
        %�ڕW�_�̌o�ܓx(phi, L), ���ʊp(alpha)�𓾂�֐�
        %Vincenty�@�̏���@ �Q�l�Fhttps://ja.wikipedia.org/wiki/Vincenty%E6%B3%95
        
        s = zeros(size(x_dist,1),1);
        for i = 1:size(x_dist,1)
            s(i,1) = norm(x_dist(i,:));
        end
        alpha_origin = pi/2 - atan2(x_dist(:,2), x_dist(:,1));
        
        sigma_origin = atan2(tan(obj.U_origin), cos(alpha_origin));
        alpha_equ = asin(cos(obj.U_origin) * sin(alpha_origin));
        u2 = cos(alpha_equ).^2 .* ((obj.a^2-obj.b^2) / obj.b^2);
        A = 1 + u2 ./ 16384 .* (4096 + u2 .* (-768 + u2 .* (320 - 175 .* u2)));
        B = u2 ./ 1024 .* (256 + u2 .* (-128 + u2 .* (74 - 47 .* u2)));
        
        sigma = s ./ (obj.b * A);
        epsilon = 1;
        i = 0;
        while epsilon >= 10^(-12)
           sigma_b = sigma;
           sigma_m = (2 * sigma_origin + sigma) / 2;
           sigma_d = B .* sin(sigma) .* cos(2 * sigma_m + 1/4 * B .*...
               (-1 + 2 * (cos(2*sigma_m)).^2 - 1/6 *B .*...
               cos(2 * sigma_m) .* (-3 + 4 * (cos(2*sigma_m)).^2)));
           sigma = s ./ (obj.b * A) + sigma_d;
           epsilon = abs((sigma - sigma_b) ./ sigma_b);
           i = i + 1;
        end
        
        phi_x = atan2(sin(obj.U_origin) .* cos(sigma) + cos(obj.U_origin) .* sin(sigma) .* cos(alpha_origin),...
            (1-obj.f) * sqrt((sin(alpha_equ)).^2 + (sin(obj.U_origin) .* sin(sigma) -...
            cos(obj.U_origin) .* cos(sigma) .* cos(alpha_origin)).^2));
        lambda = atan2(sin(sigma) .* sin(alpha_origin),...
            cos(obj.U_origin) .* cos(sigma) - sin(obj.U_origin) .* sin(sigma) .* cos(alpha_origin));
        C = obj.f / 16 * (cos(alpha_equ)).^2 .* (4 + obj.f * (4 - 3 * (cos(alpha_equ)).^2));
        L_diff = lambda - (1-C) * obj.f .* sin(alpha_equ) .*...
            (sigma + C .* sin(sigma) .* (cos(2*sigma_m) + C .*...
            cos(sigma) .* (-1 + 2 * (cos(2*sigma_m)).^2)));
        L_x = L_diff + obj.L_origin;
        alpha_x = atan2(sin(alpha_equ),...
            -sin(obj.U_origin) .* sin(sigma) + cos(obj.U_origin) .* cos(sigma) .* cos(alpha_origin));
        
        x_geo = [rad2deg(phi_x), rad2deg(L_x)];
        alpha_x = rad2deg(alpha_x);
    end
    
    function [angle_res, s] = Vincenty_inverse(obj, x_geo)
        %��_(�o�ܓx:obj.phi_origin,obj.L_origin)�ɑ΂��āA
        %�n�����W�n�C�Ӓn�_x_geo=[phi,L]�̕��ʊp(alpha)��2�_�Ԃ̋���(s)�����߂�֐�
        %alpha_origin:x_geo���猩����_�̕��ʊp, alpha:��_���猩��x_geo�̕��ʊp
        %Vincenty�@�̋t��@
        
        phi_x = deg2rad(x_geo(1));
        L_x = deg2rad(x_geo(2));
        
        if phi_x==obj.phi_origin && L_x==obj.L_origin
            angle_res = zeros(1,2);
            s = 0;
        else
            U = atan2((1-obj.f) * tan(phi_x), 1);   %�ڕW�_�̍X���ܓx
            L_diff = L_x - obj.L_origin;            %2�_�̌o�x��
            lambda = L_diff;                        %2�_�̕⏕����̌o�x��

            epsilon = 1;
            i = 0;
            while epsilon >= 10^(-12)
                lambda_b = lambda;
                s_sigma = sqrt((cos(U) * sin(lambda))^2 +...
                    (cos(obj.U_origin) * sin(U) - sin(obj.U_origin) * cos(U) * cos(lambda))^2);
                c_sigma = sin(obj.U_origin) * sin(U) + cos(obj.U_origin) * cos(U) * cos(lambda);
                sigma = atan2(s_sigma, c_sigma);
                s_alpha = cos(obj.U_origin) * cos(U) * sin(lambda) / s_sigma;
                c2_alpha = 1 - s_alpha^2;
                c_2sigma_m = c_sigma - 2 * sin(obj.U_origin) * sin(U) / c2_alpha;
                C = obj.f / 16 * c2_alpha * (4 + obj.f * (4 - 3 * c2_alpha));
                lambda = L_diff + (1 - C) * obj.f * s_alpha *...
                    (sigma + C * s_sigma * (cos(2*c_2sigma_m) +...
                    C * c_sigma * (-1 + 2 * c_2sigma_m^2)));
                epsilon = abs((lambda - lambda_b) / lambda_b);
                i = i + 1;
            end

            u2 = c2_alpha * (obj.a^2 - obj.b^2) / obj.b^2;
            A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
            B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
            sigma_d = B * s_sigma * (c_2sigma_m +...
                1/4 * B * (c_sigma * (-1 + 2 * c_2sigma_m^2) -...
                1/6 * B * c_2sigma_m * (-3 + 4 * s_sigma^2) * (-3 + 4 * c_2sigma_m^2)));
            s = obj.b * A * (sigma - sigma_d);
            alpha_origin = atan2(cos(U) * sin(lambda),...
                cos(obj.U_origin) * sin(U) - sin(obj.U_origin) * cos(U) * cos(lambda));
            alpha_x = atan2(cos(obj.U_origin) * sin(lambda),...
                -sin(obj.U_origin) * cos(U) + cos(obj.U_origin) * sin(U) * cos(lambda));

            angle_res = [rad2deg(alpha_origin), rad2deg(alpha_x)];
        end
    end
    
    function [alpha, s] = Vincenty_inverse_free(obj, x1, x2)
        %��_���W��obj.phi_origin,obj.L_origin�Ƃ��āA
        %2�_�̌o�ܓx(phi1,L1,phi2,L2)����A�e�_�ł̕��ʊp(alpha1,alpha2)��2�_�Ԃ̋���(s)�����߂�֐�
        %Vincenty_inverse���g��, ���ʎO�p�`�ɋߎ�
        
        [angle_res1, s1] = obj.Vincenty_inverse(x1);
        [angle_res2, s2] = obj.Vincenty_inverse(x2);
        alpha1 = angle_res1(1);
        alpha2 = angle_res2(1);
        s = sqrt(s1^2 + s2^2 - 2 * s1 * s2 * cosd(abs(alpha1-alpha2)));
        theta = acosd((s1^2 + s^2 - s2^2) / (2 * s1 * s));
        alpha = alpha1 - theta;
    end
    
    function x_dist = Vincenty_position(obj, x_geo)
       %��_�ɑ΂���[phi1, L]�̍��W�����߂�֐�
       %Vincenty_inverse���g��
       [angle_res, s] = obj.Vincenty_inverse(x_geo);
       alpha = 90 - angle_res(1);
       x_dist = s * [cosd(alpha), sind(alpha)];
    end
    
    function [D, alpha] = Hubeny_D(obj, phi, L)
        %2�_�̌o�ܓx���狗�������߂�֐�
        %Hubeny�̌������g�p
        
        obj.phi_origin = deg2rad(obj.phi_origin);
        obj.L_origin = deg2rad(obj.L_origin);
        phi = deg2rad(phi);
        L = deg2rad(L);
        
        phi_d = phi - obj.phi_origin;                  %2�_�̈ܓx��
        L_d = L - obj.L_origin;                        %2�_�̌o�x��
        phi_avr = mean([obj.phi_origin, phi]);         %2�_�̈ܓx�̕���
        W = sqrt(1-obj.e^2*(sin(phi_avr))^2);
        M = obj.a*(1-obj.e^2)/W^3;              %�q�ߐ��ȗ����a
        N = obj.a/W;                            %�K�ѐ��ȗ����a
        D = sqrt((phi_d*M)^2+(L_d*N*cos(phi_avr))^2);       %2�_�Ԃ̋���
        alpha = atan2(L_d*N*cos(phi_avr),phi_d*M);          %���ʊp
        
        alpha = rad2deg(alpha);
    end
    
    function [phi, L] = Hubeny_C(obj, D, alpha)
        %��_�̌o�ܓx(obj.phi_origin,obj.L_origin)�Ƃ�������̋���(D),���ʊp(alpha)����
        %�ڕW�_�̌o�ܓx(phi1,L)�����߂�֐�
        %Hubeny�̌������g�p
        
        obj.phi_origin = deg2rad(obj.phi_origin);
        obj.L_origin = deg2rad(obj.L_origin);
        alpha = deg2rad(alpha);
        
        
        
    end
    
end
end