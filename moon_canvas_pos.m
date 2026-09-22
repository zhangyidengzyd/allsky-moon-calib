function [centers, alt, az, illum] = moon_canvas_pos(files, calFile, varargin)
%MOON_CANVAS_POS 逐帧月面在"250 km 地理画布"坐标下的位置（度）
%
%   [centers, alt, az, illum] = moon_canvas_pos(files, calFile, 'Zmax', 70)
%
% 为什么需要它: 月辉是**以月面为中心的散射轮廓**。若径向归一化仍以天顶为圆心,
% 月辉就变成一个方位不对称的大梯度, 拟合不掉; 把圆心搬到月面位置, 它退化成
% 一条纯径向轮廓, 一次除掉。这一步是月夜能否看泡团的关键。
%
% files: 帧文件名 cell（或夜目录路径, 会自动列目录）
% 返回 centers: N×2, 单位度, 与 bubble_night 的 'MoonCenters' 直接对接。
%       月在地平下或算不出时间的帧返回 [0 -99]（bubble_night 会退化为天顶为圆心）。
%
% 月亮地平坐标由 moon_altaz_ref.m（Meeus 47 章, 含站心视差）给出;
% 视差 ~0.95 deg 在图上约 5 px, 不能忽略。

p = inputParser;
addParameter(p, 'Zmax', 70);
addParameter(p, 'Site', [109.133, 19.526, 0.103]);   % lon, lat, h[km]
parse(p, varargin{:});
o = p.Results;

S = load(calFile); cal = S.cal;
LAT0 = cal.obs_lat; LON0 = cal.obs_lon;
RE = 6371.0; R_OBS = RE + 0.05; R_AIR = RE + 250.0;

if ischar(files) || isstring(files)
    d = dir(fullfile(char(files), '*.PNG'));
    if isempty(d), d = dir(fullfile(char(files), '*.png')); end
    files = {d.name};
end
N = numel(files);
centers = zeros(N, 2); alt = zeros(N,1); az = zeros(N,1); illum = zeros(N,1);

for i = 1:N
    tok = regexp(files{i}, '(\d{8})(\d{6})', 'tokens');
    if isempty(tok), centers(i,:) = [0 -99]; alt(i) = -99; continue; end
    t = datetime([tok{1}{1}(1:4) '-' tok{1}{1}(5:6) '-' tok{1}{1}(7:8) ' ' ...
        tok{1}{2}(1:2) ':' tok{1}{2}(3:4) ':' tok{1}{2}(5:6)], 'TimeZone', 'UTC');
    [a, z, ~, il] = moon_altaz_ref(t, o.Site(1), o.Site(2), o.Site(3));
    alt(i) = a; az(i) = z; illum(i) = il;
    if a <= 0, centers(i,:) = [0 -99]; continue; end
    g = gamma_of_z(90 - a, R_OBS, R_AIR);
    % 沿方位 az 走地心角 g 的落点
    p0 = deg2rad(LAT0); A = deg2rad(z);
    la = asin(sin(p0)*cosd(g) + cos(p0)*sind(g)*cos(A));
    lo = deg2rad(LON0) + atan2(sind(g)*sin(A)*cos(p0), cosd(g) - sin(p0)*sin(la));
    centers(i, :) = [ (rad2deg(lo) - LON0) * cosd(LAT0), rad2deg(la) - LAT0 ];
end
end

function g = gamma_of_z(z_deg, R_OBS, R_AIR)
lo = 0; hi = acosd(R_OBS / R_AIR) * (1 - 1e-12);
for k = 1:200
    mid = 0.5 * (lo + hi);
    if atan2d(R_AIR*sind(mid), R_AIR*cosd(mid) - R_OBS) < z_deg, lo = mid; else, hi = mid; end
end
g = 0.5 * (lo + hi);
end
