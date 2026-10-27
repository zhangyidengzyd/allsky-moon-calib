function [zdeg, azdeg, r] = moon_calib_unproject(cal, x, y)
%MOON_CALIB_UNPROJECT 由像面像素坐标反解天顶角与方位角 —— 这才是做 250 km 投影时真正要用的。
%
% 输入
%   cal : moon_calib 导出的标定结构（load('calibration_params.mat') 里的 calMat）
%   x,y : 像面像素坐标（**0-based**，左上角 = (0,0)）。可以是数组。
%
% 输出
%   zdeg : 天顶角 (deg)
%   azdeg: 方位角 (deg)，北 = 0、东 = 90
%   r    : 像面半径 (px)
%
% 反演的是 r = f (z + a z^3)。因为 a 很小（|a|≈0.1）且区间固定在 z∈[0, π/2]，
% 用牛顿迭代三步就收敛到 1e-12 rad；不需要 fzero（也就不用优化工具箱）。
%
% 典型用法：把 630 nm 气辉按 250 km 高度投影到地理坐标
%   [zdeg, azdeg] = moon_calib_unproject(calMat, XX, YY);   % XX/YY 是像面网格
%   [glat, glon]  = airglow_250km(calMat.site_lat, calMat.site_lon, zdeg, azdeg);
%
% ⚠ 坐标约定：输入是 0-based；若你手上是 MATLAB 的 1-based 下标，先减 1：
%      [zdeg, azdeg] = moon_calib_unproject(calMat, colIdx - 1, rowIdx - 1);

dx = x - cal.x0;
dy = y - cal.y0;
r = hypot(dx, dy);

% 方位角：以天顶像点为极点、北 = 0 东 = 90，再扣掉图像自身的旋转角 rot
azdeg = mod(rad2deg(atan2(dx, -dy)) - cal.rot, 360);

% 牛顿反解 z：  g(z) = f (z + a z^3) - r = 0
rr = r(:);
zz = nan(size(rr));
target = rr ./ cal.f;                 % = z + a z^3
k = isfinite(target) & target > 0;
zz(k) = target(k);                    % 初值：忽略三次项
for it = 1:6
    g = zz(k) + cal.a .* zz(k).^3 - target(k);
    gp = 1 + 3 * cal.a .* zz(k).^2;
    zz(k) = zz(k) - g ./ gp;
end
zdeg = reshape(rad2deg(zz), size(r));
end
