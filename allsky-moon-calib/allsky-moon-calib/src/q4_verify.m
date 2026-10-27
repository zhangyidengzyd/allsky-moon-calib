function out = q4_verify(varargin)
%Q4_VERIFY  ④ 自检：路径体检 → 静态检查 → 投影端到端正对照 → 月亮链接通性
%
%   q4_verify
%   q4_verify('Night', '20250920', 'SkipSlow', true)
%
% 任何时候改了路径、换了定标、或怀疑"结果可疑"时跑一次。五项：
%   ① q_cfg('check')          路径存在性 + 目录分层体检
%   ② checkcode               全部脚本静态检查（要求 0 告警）
%   ③ 投影端到端正对照  ★最有判别力★
%      往**真实帧**里注入 5 个"地理坐标已知"的饱和亮斑（4 个方位角 @ z=30 与 1 个 @ z=65）：
%      先由第一性原理算出每个方向在像面上该落在哪个像素，写进去；跑完整条地理投影管线；
%      再从投影图里把标记读回来，与期望经纬度比。
%      **任何朝向 / 符号 / 翻转 / 径向模型 / rot 的错误，都会变成几 km 到上千 km 的偏差。**
%      判读：把「常量平移」（栅格原点约定，≤5 km 放过）与「位置相关残差」
%            （才是朝向/尺度类错误，≤2 km 放过）**分开判**。
%   ④ moon_calib_selftest     定标结构体字段 + 正反算往返 + 地平圈半径
%   ⑤ mc_pipeline_e2e_test    月亮标定 -> 管线 的接通性
%      （判别力来自"管线自己接受/拒绝"，不依赖我们写的公式）
%
% 参数（名值对）
%   'Night'     正对照用的夜，默认 '20250920'
%   'NMark'     标记半径(px)，默认 5
%   'GridMode'  正对照用的画布网格，默认 'geo'（与 q2_project 默认一致）
%   'SpanDeg'   正对照用的画布跨度(deg)，默认 11.98（与 q2_project 默认一致）
%   'SkipSlow'  跳过 ⑤（它要跑 20 帧完整管线），默认 false
%   'Verbose'   默认 true

p = inputParser;
addParameter(p, 'Night', '20250920');
addParameter(p, 'NMark', 5);
addParameter(p, 'GridMode', 'geo');
addParameter(p, 'SpanDeg', 11.98);
addParameter(p, 'SkipSlow', false);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
here = fileparts(mfilename('fullpath'));
addpath(here);

% 当前 rawRoot 下未必有默认夜（例如把 rawRoot 切到一期 2014 的数据根目录）。
% 找不到就自动挑"帧最多的一夜"并打印说明；否则 ③⑤ 会以"找不到夜目录"失败，
% 看起来像链路坏了，其实只是换了一批数据。
o.Night = resolve_night(cfg, o.Night);

out = struct();

% ==========================================================================
fprintf('\n########## ① 路径体检 ##########\n');
q_cfg('check');
out.cfg = q_cfg();

% ==========================================================================
acc = verify_checkcode(o.Verbose);      % 与 run_demo('verify') 同一份实现
out.checkcode = acc;

% ==========================================================================
fprintf('\n########## ③ 投影端到端正对照（注入已知地理坐标标记）##########\n');
R = proj_roundtrip(cfg, o.Night, o.NMark, o.GridMode, o.SpanDeg, o.Verbose);
out.roundtrip = R;

% ==========================================================================
fprintf('\n########## ④ 定标自检 moon_calib_selftest ##########\n');
try
    evalc('moon_calib_selftest');   % 它自己会打印; 这里只要不抛错
    fprintf('  [OK] 通过\n');
    out.selftest = true;
catch ME
    fprintf('  [失败] %s\n', ME.message);
    out.selftest = false;
end

% ==========================================================================
if ~o.SkipSlow
    fprintf('\n########## ⑤ 月亮标定 -> 管线 接通性 ##########\n');
    try
        mc_pipeline_e2e_test(fullfile(cfg.rawRoot, o.Night), 20, ...
            cfg.calNative, cfg.calFile);
        out.e2e = true;
    catch ME
        fprintf('  [失败] %s\n', ME.message);
        out.e2e = false;
    end
else
    fprintf('\n（已按要求跳过 ⑤）\n');
end

fprintf('\n########## 自检结束 ##########\n');
end

% ==========================================================================
function R = proj_roundtrip(cfg, night, nMark, gmode, spandeg, verbose)
% 注入已知地理坐标的饱和标记 -> 跑真实投影管线 -> 从投影图读回 -> 比偏差
KM = 111.195;
RE = 6371.0;
R_OBS = RE + 0.05;
R_AIR = RE + cfg.hShell;
lat0 = cfg.site(2); lon0 = cfg.site(1);
PAIRS = [0 30; 90 30; 180 30; 270 30; 135 65];   % [方位角, 天顶角]

S = load(cfg.calFile);
cal = S.cal;
x0 = double(cal.zenith(1)); y0 = double(cal.zenith(2));
rot = double(cal.rot);

srcDir = fullfile(cfg.rawRoot, night);
if exist(srcDir, 'dir') ~= 7
    error('q4_verify:noNight', '找不到夜目录: %s', srcDir);
end
fl = dir(fullfile(srcDir, cfg.pattern));
if isempty(fl), fl = dir(fullfile(srcDir, '*.PNG')); end
if numel(fl) < 3, error('q4_verify:fewFrames', '%s 帧太少', srcDir); end

work = fullfile(tempdir, 'q4_roundtrip');
if exist(work, 'dir'), rmdir(work, 's'); end
mkraw = fullfile(work, 'raw'); mkdir(mkraw);

% ---- 期望值 ----
expX = zeros(size(PAIRS, 1), 1); expY = expX;
pixX = expX; pixY = expX;
for k = 1:size(PAIRS, 1)
    az = PAIRS(k, 1); z = PAIRS(k, 2);
    r = polyval(cal.rz_poly, z);
    pixX(k) = x0 + r * sind(az + rot);
    pixY(k) = y0 - r * cosd(az + rot);
    g = solve_gamma(z, R_OBS, R_AIR);       % 注意: 返回 **度**
    gr = deg2rad(g);                        % 球面公式要用弧度
    p0 = deg2rad(lat0);
    la = asin(sin(p0) * cos(gr) + cos(p0) * sin(gr) * cosd(az));
    lo = deg2rad(lon0) + atan2(sin(gr) * sind(az) * cos(p0), ...
        cos(gr) - sin(p0) * sin(la));
    % 画布 x 坐标的约定随网格模式而变：
    %   'km'  模式把经度按 1/cos(lat0) 换算（等地面距离）
    %   'geo' 模式直接用经度差（等经纬度，与马欣论文表 3.1 的画幅一致）
    if strcmpi(gmode, 'km')
        ex = (rad2deg(lo) - lon0) * cosd(lat0);
    else
        ex = rad2deg(lo) - lon0;
    end
    ey = rad2deg(la) - lat0;
    expX(k) = ex; expY(k) = ey;
end

% ---- 注入 ----
for q = 1:3
    g = q_read_gray(fullfile(srcDir, fl(q).name));   % 按位深统一到 16 bit 量级
    if isempty(g)
        error('q4_verify:decode', '解码失败: %s', fl(q).name);
    end
    % 标记固定写 65535（满量程）：无论源帧是 8 bit 还是 16 bit，
    % 它都是画面里唯一的饱和物，与真实云/月无关。
    im = uint16(min(65535, max(0, round(g))));
    [ny, nx] = size(im);
    [XX, YY] = meshgrid(1:nx, 1:ny);
    for k = 1:size(PAIRS, 1)
        m = (XX - pixX(k)).^2 + (YY - pixY(k)).^2 <= nMark^2;
        im(m) = 65535;
    end
    imwrite(im, fullfile(mkraw, fl(q).name));
end
if verbose
    fprintf('  源夜 %s，取 3 帧注入 %d 个标记（半径 %d px）\n', night, size(PAIRS, 1), nMark);
end

% ---- 跑真实管线 ----
% ⚠ 这里必须用 'Threshold', 0 关掉自适应增益：
%   默认拉伸会把亮云也推到 255，全画布有约 3% 的像元饱和，
%   读质心时会被邻近的云拉偏（实测偏 1.5 画布像元 ≈ 3.3 km）。
%   关掉增益后天空约 23、云约 60，**标记是画面里唯一饱和的东西**，质心才干净。
geod = fullfile(work, 'geo');
airglow_geo_pipeline(mkraw, fullfile(work, 'proc'), cfg.calFile, cfg.hShell, ...
    lat0, lon0, 'BgWindow', 0, 'Margin', 30, 'ResDeg', cfg.resDeg, ...
    'GridMode', gmode, 'SpanDeg', spandeg, 'Axes', false, ...
    'GridStep', 5, 'GeoFolder', geod, 'SaveProc', false, 'Threshold', 0);

fg = dir(fullfile(geod, '*_geo.png'));
if isempty(fg)
    error('q4_verify:noGeo', '管线没有输出 *_geo.png');
end
rgb = imread(fullfile(fg(1).folder, fg(1).name));
if ismatrix(rgb), rgb = repmat(rgb, [1 1 3]); end
r16 = int16(rgb(:, :, 1)); g16 = int16(rgb(:, :, 2)); b16 = int16(rgb(:, :, 3));
sat = abs(r16 - g16) <= 4 & abs(g16 - b16) <= 4 & rgb(:, :, 1) >= 254;
[nr, nc] = size(sat);
half = (nc - 1) * cfg.resDeg / 2;

if verbose
    fprintf('  画布 %dx%d  半宽 %.4f deg  (%.2f km/px)  网格 %s，SpanDeg %s\n', ...
        nc, nr, half, cfg.resDeg * KM, gmode, ...
        tern2(isempty(spandeg), '自动', sprintf('%.2f', spandeg)));
    fprintf('\n  %3s %6s %6s | %11s %11s | %11s %11s | %8s %8s %9s\n', ...
        '#', 'az', 'z', '期望x_deg', '期望y_deg', '读出x_deg', '读出y_deg', ...
        'dx_px', 'dy_px', '偏差km');
end
errs = zeros(size(PAIRS, 1), 1);
for k = 1:size(PAIRS, 1)
    c0 = (expX(k) + half) / cfg.resDeg;      % 0-based 列
    r0 = (half - expY(k)) / cfg.resDeg;      % 0-based 行
    win = 20;
    c1 = max(1, floor(c0 - win) + 1); c2 = min(nc, ceil(c0 + win) + 1);
    r1 = max(1, floor(r0 - win) + 1); r2 = min(nr, ceil(r0 + win) + 1);
    sub = sat(r1:r2, c1:c2);
    if ~any(sub(:))
        errs(k) = NaN;
        fprintf('  %3d %6.0f %6.0f | 窗口内无饱和像元\n', k, PAIRS(k, 1), PAIRS(k, 2));
        continue
    end
    [ys, xs] = find(sub);
    cc = mean(xs) + c1 - 1;                  % 0-based 列
    rr = mean(ys) + r1 - 1;                  % 0-based 行
    xd = -half + cc * cfg.resDeg;
    yd = half - rr * cfg.resDeg;
    dd = hypot(xd - expX(k), yd - expY(k));
    errs(k) = dd;
    fprintf('  %3d %6.0f %6.0f | %11.5f %11.5f | %11.5f %11.5f | %8.2f %8.2f %9.2f\n', ...
        k, PAIRS(k, 1), PAIRS(k, 2), expX(k), expY(k), xd, yd, ...
        (xd - expX(k)) / cfg.resDeg, (yd - expY(k)) / cfg.resDeg, dd * KM);
end
ok = errs(isfinite(errs));
if isempty(ok)
    fprintf('\n  [失败] 一个标记都没读回来\n');
    R = struct('ok', false, 'errs', errs, 'geoDir', geod);
    return
end

% ---- 判读：把"常量平移"与"位置相关误差"分开 ----
% 为什么必须分开看：朝向翻转 / rot 符号 / 径向模型 / 画布 flip 这四类错误都会给出
% **随位置变化**的残差（不同方位、不同半径偏得不一样）。而整幅图的**常量平移**
% （所有标记偏同一个矢量）只可能来自栅格原点约定或插值半像元约定，它不影响泡团
% 形态，只影响绝对地理定位。本项目实测：常量偏移约 (+0.030,-0.013) deg ≈ 3.6 km，
% 位置相关残差 ≤0.3 km（在 400 km 视场内是 0.09%）。
dvx = zeros(numel(ok), 1); dvy = dvx;  %#ok<NASGU>
dxs = zeros(size(PAIRS, 1), 1); dys = dxs;
for k = 1:size(PAIRS, 1)
    if ~isfinite(errs(k)), continue; end
    c0 = (expX(k) + half) / cfg.resDeg;
    r0 = (half - expY(k)) / cfg.resDeg;
    win = 20;
    c1 = max(1, floor(c0 - win) + 1); c2 = min(nc, ceil(c0 + win) + 1);
    r1 = max(1, floor(r0 - win) + 1); r2 = min(nr, ceil(r0 + win) + 1);
    sub = sat(r1:r2, c1:c2);
    [ys, xs] = find(sub);
    dxs(k) = (mean(xs) + c1 - 1 - c0) * cfg.resDeg;
    dys(k) = ((half - (mean(ys) + r1 - 1) * cfg.resDeg) - expY(k));
end
g2 = isfinite(dxs) & isfinite(dys);
cx_off = mean(dxs(g2)); cy_off = mean(dys(g2));
resid = hypot(dxs(g2) - cx_off, dys(g2) - cy_off);
fprintf('\n  常量平移 = (%+.4f, %+.4f) deg = (%+.2f, %+.2f) km = (%.2f, %.2f) 画布像元\n', ...
    cx_off, cy_off, cx_off * KM, cy_off * KM, cx_off / cfg.resDeg, cy_off / cfg.resDeg);
fprintf('  位置相关残差: 中位 %.2f km  最大 %.2f km\n', median(resid) * KM, max(resid) * KM);
fprintf('  总偏差: 中位 %.2f km  最大 %.2f km\n', median(ok) * KM, max(ok) * KM);

rOK = (max(resid) * KM <= 2.0) && (hypot(cx_off, cy_off) * KM <= 5.0);
if rOK
    fprintf(['  [OK] 投影链几何正确（无朝向/符号/尺度错误）。\n' ...
        '       常量平移 %.2f km 属栅格原点约定，不随位置变化，\n' ...
        '       对泡团形态无影响，只影响绝对地理定位。\n'], hypot(cx_off, cy_off) * KM);
else
    fprintf('  [失败] 残差随位置变化 ⇒ 查 rot 符号 / 方向 / 径向模型 / 画布 flip\n');
end
R = struct('ok', rOK, 'errs', errs, 'geoDir', geod, ...
    'constOffKm', [cx_off * KM, cy_off * KM], ...
    'residMaxKm', max(resid) * KM, ...
    'medianKm', median(ok) * KM, 'maxKm', max(ok) * KM);
end

% ==========================================================================
function n = resolve_night(cfg, want)
% 正对照用的夜：want 不在当前 rawRoot 下时，自动挑"帧最多的一夜"。
% 理由：rawRoot 是唯一路径开关，可能被切到另一批数据（如一期 2014 的根目录）。
% 若这里硬用默认夜，③⑤ 会以"找不到夜目录"失败，看起来像链路坏了。
if exist(fullfile(cfg.rawRoot, want), 'dir') == 7
    n = want;
    return
end
d = dir(cfg.rawRoot);
d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
best = ''; bestN = 0;
for k = 1:numel(d)
    dd = fullfile(cfg.rawRoot, d(k).name);
    f = dir(fullfile(dd, cfg.pattern));
    if isempty(f), f = dir(fullfile(dd, '*.PNG')); end
    if isempty(f), f = dir(fullfile(dd, '*.png')); end
    if numel(f) > bestN
        bestN = numel(f); best = d(k).name;
    end
end
if bestN < 3
    error('q4_verify:noNight', ...
        ['rawRoot 下找不到含 >=3 帧的夜目录: %s\n' ...
         '  (期望 <rawRoot>\\<日期>\\*.PNG；请先跑 q_cfg(''check'') 体检路径)'], cfg.rawRoot);
end
fprintf('  [提示] 默认夜 "%s" 不在当前 rawRoot 下，自动改用 "%s"（%d 帧）。\n', ...
    want, best, bestN);
n = best;
end

% ==========================================================================
function g = solve_gamma(z_deg, R_OBS, R_AIR)
% 天顶角 -> 地心角(度)：视线与 H 球壳交点几何（单调，二分）
lo = 0; hi = acosd(R_OBS / R_AIR) * (1 - 1e-12);
for k = 1:200
    mid = 0.5 * (lo + hi);
    if atan2d(R_AIR * sind(mid), R_AIR * cosd(mid) - R_OBS) < z_deg
        lo = mid;
    else
        hi = mid;
    end
end
g = 0.5 * (lo + hi);
end

% ==========================================================================
function s = tern2(c, a, b)
if c, s = a; else, s = b; end
end
