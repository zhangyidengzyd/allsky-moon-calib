function [altDeg, azDeg, distKm, illum] = moon_altaz_ref(dtv, lon, lat, hgtKm)
%MOON_ALTAZ_REF 独立的月亮地平坐标参考实现（Meeus 第 47 章截断级数 + 站心视差）。
%
% 用途：作为"图像外的真值"，检验 airglow_geo_pipeline 投影出来的月亮位置对不对。
% 本文件与 moon_calib.m 内的同名算法是同一套标准公式的独立誊写，便于交叉核对：
% 用 calib_report.txt 里 [STAGE 4] 打印的 6 帧 alt/az 做自检（见 moon_calib_pipeline_check.m）。
%
% 输入
%   dtv    datetime（**UT**）标量或数组
%   lon    站点经度 (deg, 东正)
%   lat    站点纬度 (deg)
%   hgtKm  站点海拔 (km)
%
% 输出
%   altDeg, azDeg  地平高度角/方位角 (deg)，方位角 北=0 东=90
%   distKm         地心距 (km)
%   illum          月面被照亮比例 [0,1]

altDeg = zeros(size(dtv)); azDeg = altDeg; distKm = altDeg; illum = altDeg;
for n = 1:numel(dtv)
    jd = local_jd(dtv(n));
    [lam, bet, eps, dist] = local_moon_ecl(jd);
    [ra, de] = local_ecl2eq(lam, bet, eps);
    sinPi = 6378.14 / dist;

    phi = deg2rad(lat);
    u   = atan(0.99664719 * tan(phi));
    rs  = 0.99664719 * sin(u) + (hgtKm / 6378.14) * sin(phi);
    rc  = cos(u) + (hgtKm / 6378.14) * cos(phi);

    gst = local_gmst_deg(jd);
    H  = deg2rad(gst + lon) - ra;
    dA = atan2(-rc * sinPi * sin(H), cos(de) - rc * sinPi * cos(H));
    raT = ra + dA;
    deT = atan2((sin(de) - rs * sinPi) * cos(dA), cos(de) - rc * sinPi * cos(H));

    Ht = deg2rad(gst + lon) - raT;
    z  = acos(max(-1, min(1, sin(phi) * sin(deT) + cos(phi) * cos(deT) * cos(Ht))));
    A  = atan2(-sin(Ht) * cos(deT), sin(deT) * cos(phi) - cos(deT) * sin(phi) * cos(Ht));

    altDeg(n) = 90 - rad2deg(z);
    azDeg(n)  = mod(rad2deg(A), 360);
    distKm(n) = dist;
    illum(n)  = local_illum(lam, bet, jd);
end
end

% ---------------------------------------------------------------------------
function jd = local_jd(t)
y = year(t); m = month(t);
d = day(t) + (hour(t) + minute(t)/60 + second(t)/3600)/24;
if m <= 2, y = y - 1; m = m + 12; end
A = floor(y/100); B = 2 - A + floor(A/4);
jd = floor(365.25*(y+4716)) + floor(30.6001*(m+1)) + d + B - 1524.5;
end

function g = local_gmst_deg(jd)
T = (jd - 2451545.0)/36525.0;
g = mod(280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*T*T, 360.0);
end

function [lam, bet, eps, distKm] = local_moon_ecl(jd)
T  = (jd - 2451545.0)/36525.0;
r  = @deg2rad;
Lp = r(mod(218.3164477 + 481267.88123421 * T, 360));
D  = r(mod(297.8501921 + 445267.1114034  * T, 360));
M  = r(mod(357.5291092 + 35999.0502909  * T, 360));
Mp = r(mod(134.9633964 + 477198.8675055 * T, 360));
F  = r(mod(93.2720950  + 483202.0175233 * T, 360));
E  = 1.0 - 0.002516*T - 0.0000074*T*T;

lam = Lp + r(6.288774)*sin(Mp) + r(1.274027)*sin(2*D - Mp) ...
        + r(0.658314)*sin(2*D) + r(0.213618)*sin(2*Mp) ...
        - r(0.185116)*E*sin(M) - r(0.114332)*sin(2*F) ...
        + r(0.058793)*sin(2*D - 2*Mp) + r(0.057066)*E*sin(2*D - M - Mp) ...
        + r(0.053322)*sin(2*D + Mp) + r(0.045758)*E*sin(2*D - M) ...
        - r(0.040923)*E*sin(M - Mp) - r(0.034720)*sin(D) ...
        - r(0.030383)*E*sin(M + Mp) + r(0.015327)*sin(2*D - 2*F) ...
        - r(0.012528)*sin(Mp + 2*F) + r(0.010980)*sin(Mp - 2*F);

bet = r(5.128622)*sin(F) + r(0.280602)*sin(Mp + F) + r(0.277693)*sin(Mp - F) ...
        + r(0.173237)*sin(2*D - F) + r(0.055413)*sin(2*D - Mp + F) ...
        + r(0.046271)*sin(2*D - Mp - F) + r(0.032573)*sin(2*D + F) ...
        + r(0.017198)*sin(2*Mp + F);

distKm = 385000.56 - 20905.355*cos(Mp) - 3699.111*cos(2*D - Mp) ...
         - 2955.968*cos(2*D) - 569.925*cos(2*Mp);

A1 = r(mod(119.75 + 131.849   * T, 360));
A2 = r(mod(53.09  + 479264.290 * T, 360));
A3 = r(mod(313.45 + 481266.484 * T, 360));
lam = lam + r(0.004*(sin(A1) + sin(Lp - F)) + 0.000318*sin(A2));
bet = bet + r(0.002*(sin(A3) + sin(A1 - F) + sin(A1 + F) + sin(Lp - Mp) - sin(Lp + Mp)));

eps = r(23.4392911 - 0.0130042*T);
end

function [ra, de] = local_ecl2eq(lam, bet, eps)
ra = atan2(sin(lam)*cos(eps) - tan(bet)*sin(eps), cos(lam));
de = asin(sin(bet)*cos(eps) + cos(bet)*sin(eps)*sin(lam));
ra = mod(ra, 2*pi);
end

function k = local_illum(lamM, betM, jd)
T = (jd - 2451545.0)/36525.0;
L0 = 280.46646 + 36000.76983*T + 0.0003032*T*T;
M  = deg2rad(mod(357.52911 + 35999.05029*T - 0.0001537*T*T, 360));
C  = (1.914602 - 0.004817*T - 0.000014*T*T)*sin(M) ...
   + (0.019993 - 0.000101*T)*sin(2*M) + 0.000289*sin(3*M);
lamS = deg2rad(mod(L0 + C, 360.0));
k = 0.5*(1.0 + (-cos(lamM - lamS)*cos(betM)));
end
