classdef lon_lat
   %方位角・距離⇔経緯度間の変換クラス
   %地球楕円体はWGS84(GooglマップやGPSと同じ)を採用
   %方位角は 北0deg, CW (推測)
   
properties
    a = 6378137.06;         %地球赤道半径 [m] (WGS84)
    b = 0;                  %地球極半径 [m]
    f = 1/298.257223563;    %地球の扁平率 (WGS84)
    e = 0;                  %楕円断面の離心率
    phi_origin = 0;         %基準点緯度 [deg]
    L_origin = 0;           %基準点経度 [deg]
    U_origin = 0;           %基準点の更成緯度 [deg]
end

methods
    function obj = lon_lat(pos)     %コンストラクタメソッド
        obj.b = obj.a*(1-obj.f);
        obj.e = sqrt(obj.f*(2-obj.f));
        obj.phi_origin = deg2rad(pos(1));
        obj.L_origin = deg2rad(pos(2));
        obj.U_origin = atan2((1-obj.f) * tan(obj.phi_origin),1);
    end
    
    function [x_geo, alpha_x] = Vincenty_direct(obj, x_dist)
        %基準点の経緯度(obj.phi_origin, obj.L_origin), 目標点座標x_dist = [x,y]から
        %目標点の経緯度(phi, L), 方位角(alpha)を得る関数
        %Vincenty法の順解法 参考：https://ja.wikipedia.org/wiki/Vincenty%E6%B3%95
        
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
        %基準点(経緯度:obj.phi_origin,obj.L_origin)に対して、
        %地理座標系任意地点x_geo=[phi,L]の方位角(alpha)と2点間の距離(s)を求める関数
        %alpha_origin:x_geoから見た基準点の方位角, alpha:基準点から見たx_geoの方位角
        %Vincenty法の逆解法
        
        phi_x = deg2rad(x_geo(1));
        L_x = deg2rad(x_geo(2));
        
        if phi_x==obj.phi_origin && L_x==obj.L_origin
            angle_res = zeros(1,2);
            s = 0;
        else
            U = atan2((1-obj.f) * tan(phi_x), 1);   %目標点の更成緯度
            L_diff = L_x - obj.L_origin;            %2点の経度差
            lambda = L_diff;                        %2点の補助球上の経度差

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
        %基準点座標をobj.phi_origin,obj.L_originとして、
        %2点の経緯度(phi1,L1,phi2,L2)から、各点での方位角(alpha1,alpha2)と2点間の距離(s)を求める関数
        %Vincenty_inverseを拡張, 平面三角形に近似
        
        [angle_res1, s1] = obj.Vincenty_inverse(x1);
        [angle_res2, s2] = obj.Vincenty_inverse(x2);
        alpha1 = angle_res1(1);
        alpha2 = angle_res2(1);
        s = sqrt(s1^2 + s2^2 - 2 * s1 * s2 * cosd(abs(alpha1-alpha2)));
        theta = acosd((s1^2 + s^2 - s2^2) / (2 * s1 * s));
        alpha = alpha1 - theta;
    end
    
    function x_dist = Vincenty_position(obj, x_geo)
       %基準点に対する[phi1, L]の座標を求める関数
       %Vincenty_inverseを拡張
       [angle_res, s] = obj.Vincenty_inverse(x_geo);
       alpha = 90 - angle_res(1);
       x_dist = s * [cosd(alpha), sind(alpha)];
    end
    
    function [D, alpha] = Hubeny_D(obj, phi, L)
        %2点の経緯度から距離を求める関数
        %Hubenyの公式を使用
        
        obj.phi_origin = deg2rad(obj.phi_origin);
        obj.L_origin = deg2rad(obj.L_origin);
        phi = deg2rad(phi);
        L = deg2rad(L);
        
        phi_d = phi - obj.phi_origin;                  %2点の緯度差
        L_d = L - obj.L_origin;                        %2点の経度差
        phi_avr = mean([obj.phi_origin, phi]);         %2点の緯度の平均
        W = sqrt(1-obj.e^2*(sin(phi_avr))^2);
        M = obj.a*(1-obj.e^2)/W^3;              %子午線曲率半径
        N = obj.a/W;                            %卯酉線曲率半径
        D = sqrt((phi_d*M)^2+(L_d*N*cos(phi_avr))^2);       %2点間の距離
        alpha = atan2(L_d*N*cos(phi_avr),phi_d*M);          %方位角
        
        alpha = rad2deg(alpha);
    end
    
    function [phi, L] = Hubeny_C(obj, D, alpha)
        %基準点の経緯度(obj.phi_origin,obj.L_origin)とそこからの距離(D),方位角(alpha)から
        %目標点の経緯度(phi1,L)を求める関数
        %Hubenyの公式を使用
        
        obj.phi_origin = deg2rad(obj.phi_origin);
        obj.L_origin = deg2rad(obj.L_origin);
        alpha = deg2rad(alpha);
        
        
        
    end
    
end
end