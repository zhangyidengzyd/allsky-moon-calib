function [altDeg, azDeg] = sun_altaz_ref(dtv, lon, lat)
%SUN_ALTAZ_REF 独立的太阳地平坐标参考实现（Meeus 第 25 章低精度 + 黄赤转换）。
%
%   [alt, az] = sun_altaz_ref(datetime(2014,9,30,12,0,0), 109.133, 19.526)
%
% 用途
% ----
% 与 moon_altaz_ref 对称：给出**图像外的真值**，供两条用途：
%   1) 判定一帧是不是"暮光帧"。8 bit 的一期 2014 数据里，太阳只要高于约 -12 deg，
%      大片天光就能顶到满量程，被 moon_calib 的 [STAGE 2] 当成"饱和月核"。
%      这些假目标的质心跟着**太阳**走而不是月亮 ⇒ 必须按太阳高度把它们剔掉。
%      （moon_calib 内部的 mc_sun_altaz 是同一套公式的另一份誊写，互为对照。）
%   2) 事后核对：calib_2014 用本函数，moon_calib 用 mc_sun_altaz，两者的差应当 <0.05 deg。
%
% 输入
%   dtv    datetime（**UT**）标量或数组
%   lon    站点经度 (deg, 东正)
%   lat    站点纬度 (deg)
%
% 输出
%   altDeg, azDeg  地平高度角/方位角 (deg)，方位角 北=0 东=90（与像面 phi 约定一致）
%
% 精度：太阳黄经 ~0.01 deg（Meeus 低精度级数），对应高度角误差 <0.02 deg —— 远小于
%       本用途需要的 1 deg 量级判据。太阳视差 ~8.8" 忽略，不做站心修正。

n = numel(dtv);
altDeg = zeros(size(dtv));
azDeg = altDeg;
for k = 1:n
    jd = local_jd(dtv(k));
    T = (jd - 2451545.0) / 36525.0;

    % --- 太阳几何黄经（Meeus 25.2-25.4）---
    L0 = mod(280.46646 + 36000.76983 * T + 0.0003032 * T * T, 360);
    M = deg2rad(mod(357.52911 + 35999.05029 * T - 0.0001537 * T * T, 360));
    C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * sin(M) ...
      + (0.019993 - 0.000101 * T) * sin(2 * M) ...
      + 0.000289 * sin(3 * M);
    lam = deg2rad(mod(L0 + C, 360));
    eps = deg2rad(23.4392911 - 0.0130042 * T);

    % --- 黄道 -> 赤道 ---
    ra = atan2(sin(lam) * cos(eps), cos(lam));
    de = asin(sin(lam) * sin(eps));

    % --- 赤道 -> 地平（无站心视差）---
    gst = mod(280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T * T, 360);
    H = deg2rad(gst + lon) - ra;
    phi = deg2rad(lat);
    z = acos(max(-1.0, min(1.0, sin(phi) * sin(de) + cos(phi) * cos(de) * cos(H))));
    A = atan2(-sin(H) * cos(de), sin(de) * cos(phi) - cos(de) * sin(phi) * cos(H));

    altDeg(k) = 90.0 - rad2deg(z);
    azDeg(k) = mod(rad2deg(A), 360.0);
end
end

% ---------------------------------------------------------------------------
function jd = local_jd(t)
y = year(t); m = month(t);
d = day(t) + (hour(t) + minute(t) / 60 + second(t) / 3600) / 24;
if m <= 2, y = y - 1; m = m + 12; end
A = floor(y / 100); B = 2 - A + floor(A / 4);
jd = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + B - 1524.5;
end
