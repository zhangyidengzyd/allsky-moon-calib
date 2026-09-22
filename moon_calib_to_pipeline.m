function cal = moon_calib_to_pipeline(calMat, outFile, varargin)
%MOON_CALIB_TO_PIPELINE 把 moon_calib 的标定转成 airglow_geo_pipeline 能直接 load 的 cal 结构。
%
%   cal = moon_calib_to_pipeline(calMat)                       % 只返回结构
%   cal = moon_calib_to_pipeline(calMat, outFile)              % 同时 save(outFile,'cal')
%   cal = moon_calib_to_pipeline(calMat, outFile, 'Radius', 512.93)
%
% 参数
%   calMat   moon_calib 的标定结构（load('calibration_params.mat') 的 calMat）
%   outFile  输出 .mat 路径（可省略/给 ''）。写入变量名固定为 `cal`
%   'Radius'   写进 cal.radius 的视场边缘半径 (px)，默认 512.93（实测视场软边缘半高点）
%   'Zenith'   手工覆盖天顶像点（1-based px）；默认由 calMat 的 0-based 值 +1
%   'Verbose'  打印换算明细，默认 true
%
% ---------------------------------------------------------------------------
% 为什么需要这个转换（每条都是实测确认过的约定，改错任何一条都是系统性偏差）
% ---------------------------------------------------------------------------
% (1) 坐标基准：moon_calib 内部全程 **0-based**（左上角像元 = (0,0)）；
%     而 airglow_geo_pipeline 里写的是
%         [X, Y] = meshgrid(1:w, 1:h);  mask = hypot(X-cx, Y-cy) <= r_use;
%     即 cx/cy 是 **1-based**。所以 zenith 必须 +1。
%     不 +1 的话天顶会整体错 1 px（0.14 deg），小但白错。
%
% (2) rot 的符号：**不用翻**。moon_calib 的 rot 定义是
%         x = x0 + r*sin(A + rot),   y = y0 - r*cos(A + rot)
%     用户 run.m 管线里读进来的那个文件自己写着（字段 rot_convention）：
%         'star-fit parameter in x=zx+r*sin(A+rot), y=zy-r*cos(A+rot);
%          downstream affine2d correction uses theta=+rot'
%     两者定义完全一致 ⇒ cal.rot = calMat.rot 原样传入。
%     （推导一下更放心：管线算 X_dest = c + r(sin A, -cos A)，再
%       X_src = transformPointsInverse(tform_rot, X_dest)。affine2d 用行向量约定，
%       T2 = [cos t, -sin t; sin t, cos t] 作用在 (u,v) 上得到
%       (u cos t + v sin t, -u sin t + v cos t) = R(-t) 作用在列向量上。
%       于是 X_src = c + R(+t)·r(sin A, -cos A) = c + r(sin(A+t), -cos(A+t))，
%       即从源图采到的方位是 A+t。源图里方位 A' 的像点在 A'+rot 处
%       ⇒ 采到的天空方位 = (A+t) - rot。要它等于 A 就必须 t = +rot。✔）
%
% (3) 径向多项式单位：管线里 make_rz_eval 用
%         if abs(rz(end-1)) > 50 -> 系数按 **弧度**，polyval(rz, deg2rad(zdeg))
%         else                  -> 系数按 **度**，  polyval(rz, zdeg)
%     rz(end-1) 永远是 z^1 的系数（末位是 z^0）。本项目 r = f(z+a z^3) 若用度表示，
%     z^1 系数 = f*pi/180 = 7.37 < 50 ⇒ 走"度"分支。所以这里输出 **度** 系数、8 元
%     （和用户现有 calibration_params.mat 的格式一致：rz_poly_units = 'degree'）：
%         r_px = c5*zdeg^3 + c1*zdeg,  c1 = f*pi/180,  c5 = f*a*(pi/180)^3
%     验算：zdeg=90 -> 7.3696*90 - 2.2416e-4*729000 = 663.26 - 163.41 = 499.85 px ✔
%     （若改用弧度系数也行：z^1 系数 = f = 422 > 50 ⇒ 自动走弧度分支；
%       但为和现有文件同格式，这里默认给"度"。）
%
% (4) obs_lat/obs_lon：管线里 load_calibration 之后有
%         if isfield(cal,'obs_lat') && ... -> obs_lat = cal.obs_lat; obs_lon = cal.obs_lon;
%     也就是说 **cal 里的站点坐标会覆盖 run.m 传进来的实参**。
%     本项目官方站点是 (109.133E, 19.526N)，而 run.m 现在硬写 (19.5, 109.2)。
%     本函数默认写官方值（差 0.067 deg 经度 ⇒ 250 km 上约 7 km），
%     如果要保留 run.m 的入参行为，用 'SiteMode','keep'。
%
% (5) 手性：本数据实测为 **直接（无镜像）**（rms 2.11 px vs 镜像 23.4 px）。
%     管线没有镜像开关，直接手性是它预期的约定，因此无需额外处理。
%
% 用法（最短路径）
%   S = load('mooncal_work\mooncal_out_matlab\calibration_params.mat');
%   cal = moon_calib_to_pipeline(S.calMat, 'E:\qihui2\Processed_Output\calibration_params_moon.mat');
%   % 然后把 run.m 第 3 个实参换成这个新文件即可。

p = inputParser;
addParameter(p, 'Radius',   512.93, @(v) isscalar(v) && v > 0);
addParameter(p, 'Zenith',   [],     @(v) isempty(v) || numel(v) == 2);
addParameter(p, 'SiteMode', 'official', @(v) ischar(v) && any(strcmpi(v, {'official', 'keep'})));
addParameter(p, 'Verbose',  true,   @(v) islogical(v) || isnumeric(v));
parse(p, varargin{:});
o = p.Results;

%% ---------- 1) 天顶像点：0-based -> 1-based ----------
if isempty(o.Zenith)
    cal.zenith = [calMat.x0 + 1, calMat.y0 + 1];
else
    cal.zenith = o.Zenith(:).';
end

%% ---------- 2) rot：直接沿用（定义一致）----------
cal.rot = calMat.rot;

%% ---------- 3) 径向多项式：换成“度”系数的 8 元向量 ----------
DEG = pi / 180;
f = calMat.f;
a = calMat.a;
c1 = f * DEG;                 % z^1
c5 = f * a * DEG^3;           % z^3

% polyval 降幂排列： z^7 z^6 z^5 z^4 z^3 z^2 z^1 z^0
cal.rz_poly = [0, 0, 0, 0, c5, 0, c1, 0];

% 自检：多项式必须能重现 f(z + a z^3)
zt = [0.002, 0.5, 1.0, pi/2];
r_src = f * (zt + a * zt.^3);
r_dst = polyval(cal.rz_poly, rad2deg(zt));
assert(max(abs(r_src - r_dst)) < 1e-6, ...
    'moon_calib:polyMismatch', 'rz_poly 转换自检失败（最大差 %.3g px）', max(abs(r_src - r_dst)));
% 断言启发式分支确实落在“度”分支上
assert(abs(cal.rz_poly(end-1)) <= 50, ...
    'moon_calib:unitBranch', 'z^1 系数 %.4g 会触发管线的“弧度”分支，单位约定不一致', cal.rz_poly(end-1));

%% ---------- 4) radius / 站点坐标 ----------
cal.radius = o.Radius;
% ⚠ 这里必须给个默认值：下面第 8 步（附带信息）会无条件读 cal.rim_px。
%   历史 bug：只有 calMat.rim_px 有限时才赋值，于是用「自己合成的 calMat」
%   （例如 calib_2014 的"纯方位圆心 + 固定圆心"采用解）调用本函数时，
%   会在 cal.R_measured = cal.rim_px 处报「无法识别的字段 "rim_px"」。
cal.rim_px = NaN;
if isfield(calMat, 'rim_px') && ~isempty(calMat.rim_px) && isfinite(calMat.rim_px)
    cal.rim_px = calMat.rim_px;
end

switch lower(o.SiteMode)
    case 'official'
        cal.obs_lat = 19.526;
        cal.obs_lon = 109.133;
    case 'keep'
        if isfield(calMat, 'site_lat')
            cal.obs_lat = calMat.site_lat;
            cal.obs_lon = calMat.site_lon;
        else
            error('moon_calib:noSite', 'calMat 里没有 site_lat/site_lon，无法 SiteMode=keep');
        end
end

%% ---------- 5) 附带信息（不参与计算，便于日后核对来源）----------
cal.r90            = calMat.f * (pi/2 + calMat.a * (pi/2)^3);
cal.R_measured     = cal.rim_px;
% ★ 2026-09-17 新增：盘边对应的天顶角。**它不是 90 deg** —— 仪器标称 180 deg 视场、
%   有效约 160 deg；一期实测盘边 ≈87 deg（实测锚点只到 z=76 deg）。因此
%   "f = R_rim/(pi/2)" 这类换算不成立（一期会低估 f 约 10%）。见 calib_2014.m (E5) / fov_check.m。
cal.z_rim_deg      = NaN;
if isfield(calMat, 'z_rim_deg') && ~isempty(calMat.z_rim_deg) && isfinite(calMat.z_rim_deg)
    cal.z_rim_deg = calMat.z_rim_deg;
end
cal.fov_note       = '';
if isfield(calMat, 'fov_note') && ~isempty(calMat.fov_note)
    cal.fov_note = calMat.fov_note;
end
cal.rz_poly_units  = 'degree';
cal.rz_poly_model  = 'r_px = polyval(rz_poly, zenith_angle_deg)';
cal.rot_units      = 'degree, signed (-180,180]';
cal.rot_convention = ['x=zx+r*sin(A+rot), y=zy-r*cos(A+rot); ' ...
                      'downstream affine2d correction uses theta=+rot'];
cal.azimuth_offset = 0;
cal.model          = sprintf('equid+ (moon calib): r = %.6g * (z + %.6g z^3), z in rad', f, a);
cal.focal          = f;
cal.a              = a;
cal.rms_px         = calMat.rms_px;
cal.n_frames       = calMat.n_frames;
cal.n_nights       = calMat.n_nights;
cal.source         = 'moon_calib (moon-based orientation + radial calibration)';
cal.version        = 'moon_calib_m1 -> airglow_geo_pipeline adapter';

%% ---------- 6) 健全性断言（这几条任意一条不过就不该交给投影管线）----------
assert(isscalar(cal.rot) && isfinite(cal.rot), 'rot 非有限标量');
assert(all(isfinite(cal.zenith)) && all(cal.zenith >= 1), 'zenith 非法（必须 1-based 且为正）');
assert(numel(cal.rz_poly) == 8 && all(isfinite(cal.rz_poly)), 'rz_poly 必须 8 元有限值');
r_rim_test = polyval(cal.rz_poly, 90);
assert(r_rim_test > 400 && r_rim_test < 600, 'r(90 deg) = %.1f px 落在合理区间之外', r_rim_test);

%% ---------- 7) 写盘 ----------
if nargin >= 2 && ~isempty(outFile)
    save(outFile, 'cal');
    if o.Verbose
        fprintf('已写 -> %s\n', outFile);
    end
end

if o.Verbose
    fprintf('\n=== moon_calib -> airglow_geo_pipeline 适配结果 ===\n');
    fprintf('  zenith  (1-based)  = (%.3f, %.3f)   [calMat 0-based: (%.3f, %.3f) +1]\n', ...
        cal.zenith(1), cal.zenith(2), calMat.x0, calMat.y0);
    fprintf('  rot                = %+.4f deg        [原样沿用，theta = +rot]\n', cal.rot);
    fprintf('  radius             = %.2f px\n', cal.radius);
    fprintf('  rz_poly (deg, 8)   = [%s]\n', strrep(sprintf('%.6g ', cal.rz_poly), '  ', ' '));
    fprintf('  r(90 deg)          = %.2f px  (= r_horizon)\n', cal.r90);
    fprintf('  obs_lat/obs_lon    = %.3f, %.3f   [SiteMode = %s]\n', ...
        cal.obs_lat, cal.obs_lon, o.SiteMode);
    fprintf('  多项式自检          : 与 f(z+a z^3) 最大差 %.2e px\n', max(abs(r_src - r_dst)));
    fprintf('==================================================\n\n');
end
end
