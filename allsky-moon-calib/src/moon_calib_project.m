function [x, y, zdeg, r] = moon_calib_project(cal, lat, dec, Hdeg, azdeg)
%MOON_CALIB_PROJECT 用 moon_calib 解出的标定，把天球坐标正向投影到像面像素坐标。
%
% 输入（除 cal 外都可以是同样大小的数组）
%   cal  : moon_calib 导出的标定结构（load('calibration_params.mat') 里的 calMat）
%   lat  : 站点地理纬度 (deg)
%   dec  : 目标赤纬 (deg)
%   Hdeg : 目标时角 (deg)，H = LST - RA
%   azdeg: (可选) 目标方位角 (deg)，北 = 0、东 = 90。
%          不给就由 (lat, dec, H) 现算 —— 给的方位角必须与它一致（差一个常数 rot 会被
%          直接吸收，所以一般建议不给，让公式统一算出，避免两处约定不一致）。
%
% 输出（**0-based 像素坐标**，左上角像元 = (0,0)、右下角 = (1023,1023)）
%   x, y : 像面像素坐标
%   zdeg : 天顶角 (deg)
%   r    : 像面半径 (px)
%
% 用法
%   calMat = load('calibration_params.mat'); calMat = calMat.calMat;
%   [az, ~] = moon_calib_az(35.0, dec, H);       % 只是演示，实际见下
%   [x, y] = moon_calib_project(calMat, 19.526, dec, H);   % az 省略，自动算
%
% ⚠ 三个必须记住的约定（错了就是系统性偏差）
%   1) 输出是 0-based。要拿去索引 MATLAB 数组（1-based）必须 +1：
%        I(round(y)+1, round(x)+1)
%      这是本项目最容易踩的一条 —— 标定值 x0≈538.9 / y0≈476.5 都是 0-based。
%   2) 方位角北 = 0、东 = 90（不是数学上的东 = 0 逆时针）。
%   3) z 的单位，公式里用**弧度**、报告里用度。函数内部已经处理好了。

if nargin < 5 || isempty(azdeg)
    azdeg = moon_calib_az(lat, dec, Hdeg);
end

% 天顶角：cos z = sin(lat)sin(dec) + cos(lat)cos(dec)cos(H)
cz = sind(lat) .* sind(dec) + cosd(lat) .* cosd(dec) .* cosd(Hdeg);
cz = min(1, max(-1, cz));
zz = acos(cz);                       % rad

% 径向投影  r = f (z + a z^3)   （本项目最佳模型 'equid+'）
r = cal.f .* (zz + cal.a .* zz.^3);

% 像面坐标（0-based）
x = cal.x0 + r .* sin(deg2rad(azdeg) + deg2rad(cal.rot));
y = cal.y0 - r .* cos(deg2rad(azdeg) + deg2rad(cal.rot));

zdeg = rad2deg(zz);
end


function azdeg = moon_calib_az(lat, dec, Hdeg)
%MOON_CALIB_AZ 由 (lat, dec, H) 算方位角 (deg)：北 = 0、东 = 90。
num = -sind(Hdeg) .* cosd(dec);
den = sind(dec) .* cosd(lat) - cosd(dec) .* sind(lat) .* cosd(Hdeg);
azdeg = mod(rad2deg(atan2(num, den)), 360);
end
