function zen = zenith_from_rim(nightDir, varargin)
%ZENITH_FROM_RIM  Boresight point + field-edge radius from the circular field rim.
%
%   本函数**不需要星**。旋转对称光学系统里 r 只是"离轴角"的函数, 所以 r = const 的
%   轨线是**以光轴与探测器交点(boresight 像素)为圆心的圆**; 视场边界(天光跌到盲电平
%   的半径处)就是这样一个圆 ⇒ 圆心 = boresight 像素。
%
%   ★★★ 先读这一节, 它决定这个输出能怎么用 ★★★
%
%   [结论 1] 圆心是 **boresight**, 不是"天顶像点"。
%       圆心 = 光轴 ∩ 探测器。它等于天顶像点的**前提**是: 该时刻光轴指向天顶。
%       本仪器 = 儋州站双通道全天空气辉成像仪 DCAI, **固定式、无扫描机构**,
%       光轴相对天顶的姿态完全由机械安装/调平决定, **没有"驻停位置"这种可变指向状态**;
%       指向是否朝天顶是一个**安装事实**, 而**不能靠图像自证**:
%       光轴偏天顶 t 对半径的一阶效应是方位一次谐波
%           dr = -f (1 + 3 a z^2) * t * cos(phi - phi_t)
%       而圆心平移产生的是**同频率、振幅不随 z 变**的谐波 ⇒ 两者只在"振幅随 z 的
%       变化"上可分。样本方位覆盖 <30 deg 时, 倾斜几乎完全被 (圆心, f) 吸收。
%       实测(本项目 349 帧月亮样本, 方位只覆盖 24 deg): 把倾斜 t 作自由参数扫一遍,
%       rms 从 2.120 px (t=0) 单调降到 1.742 px (t=15 deg), 同时圆心滑走 100+ px、
%       f 变化 4% —— 这是**平谷退化**, 不是"测到了倾斜"。留一夜交叉给出
%       12/12/12/12/2 deg 的不稳定答案。
%       ⇒ 唯一的硬约束来自外部: 本函数测到的圆心就是 boresight, 拿它去和
%         "靠月亮/星定出的天顶像点"比较, 差值 = f*t, 才给出 t 的真上限。
%         用 'ZenithRef' 传入那个天顶像点即可自动给出判定。
%       实测: 边界圆圆心 (543.82, 474.35) vs 月亮定标圆心 (538.89, 476.47),
%       差 5.37 px = 0.73 deg (以 f = 422.25 px/rad 折算) ⇒ 这批观测的 boresight
%       确实基本朝天顶, **但这是"测出来"的, 不是"假设"的**; 换一批数据、或相机重新
%       安装/调平过, 都必须重测。
%       (注: "扫描镜方位 90 deg / 仰角 45 deg" 是同站 DCOI 双通道光学干涉仪的参数,
%        与本全天相机无关 —— 2026-09-16 用户更正。)
%
%   [结论 2] 半径 R_rim 是**无假设的测量**; 但 "f = R/theta" **不是测量**。
%       视场边界对应的天顶角 theta_rim 由光学设计决定, 可能是 77 deg 也可能是 90 deg,
%       图像本身无法区分。f = R/theta 对 theta 是 1/theta 敏感:
%           theta_rim = 75 deg -> f 比 theta=90 deg 时大 20%
%       ⇒ **不要**用 FovHalfDeg 默认为 90 的 focal 去喂标定链。本函数已把
%         FovHalfDeg 默认值改为 NaN: 不给就只输出圆心与半径, focal = NaN,
%         cal.rz_poly 全 NaN。若确实要写入一个"按假设换算"的尺度, 必须显式给
%         'FovHalfDeg' 且 'ExportScale', true, 且结果会带上
%         cal.focalIsMeasurement = false 与来源字符串 —— 它**不是**测量值。
%       这正是本项目 v6 脚本 [发现B] 的结论(当时把 r(90)=R_measured 这条硬锚定
%       改成可选且默认关闭): 若可用视场只有 ~155 deg, 亮圈对应 z≈77 deg,
%       那条约束会把解空间拖到一个"妥协解"上。径向尺度的正确来源是月亮/星定标。
%
%   USAGE
%     zen = zenith_from_rim(nightDir)
%     zen = zenith_from_rim(nightDir, 'ZenithRef', [538.89 476.47])
%     zen = zenith_from_rim(nightDir, 'FovHalfDeg', 90, 'ExportScale', true)
%
%   INPUT
%     nightDir : folder (or cell array of folders) holding the L0 PNG frames.
%
%   OPTIONS (name/value)
%     'DarkFile'      dark/bias frame path, subtracted before fitting ('' = none)
%     'NFrames'       frames stacked into the median, default 30 (max 60)
%     'FovHalfDeg'    zenith angle AT the field edge, degrees. **default NaN = 未知**.
%                     只在明确知道光学设计时给值; 90 => 假设 180 deg 视场。
%     'HalfFovRange'  [lo hi] deg, 只用于打印"theta 未知时 f 的取值范围"。默认 [70 110]
%     'ExportScale'   logical, default false. 只有 true 且 FovHalfDeg 有限时,
%                     focal 与 rz_poly 才会写进 cal; 否则写 NaN。
%     'ZenithRef'     已知的天顶像点 [x y] (1-based, 与 zen.center 同约定)。
%                     给了就做"圆心=天顶"前提检验并打印 / 判定。
%     'ScalePxRad'    用来把像素差折成角度的径向尺度 (px/rad)。**应当传月亮/星定标值**
%                     (本项目 422.25)。不传则临时用 R_rim/90deg 折算, 并标注为粗糙。
%     'SearchHalf'    half-width of the half-maximum search band, px, default 100
%     'NSig'          sigma clipping in the circle fit, default 2.5
%     'SaveMat'       path of the output .mat ('' = skip)
%     'SavePreview'   path of the preview png ('' = skip)
%
%   OUTPUT struct `zen`
%     center        [x0 y0] in 1-based MATLAB pixel coords == **boresight** 像素
%     center0       same, in 0-based (x,y) like OpenCV / python
%     radius        field-edge radius, px            <-- 无假设测量
%     rms           circle-fit residual, px
%     nPoints       edge samples retained
%     focal         px/rad. **NaN 除非显式 'FovHalfDeg' + 'ExportScale'**。
%                   即便给了值, focalIsMeasurement 仍为 false。
%     centerHalf    centre from an independent half-maximum edge detector
%     dRadius       radius difference between the two detectors (systematic check)
%     zenithCheck   仅在给 'ZenithRef' 时存在: 结构体, 含 d_px / d_deg / verdict
%     cal           struct ready to be merged into calibration_params.mat
%     rotKnown      false -- this method CANNOT determine the rotation angle
%
%   IMPORTANT LIMITS
%     * 圆心 = boresight: 确定。圆心 = 天顶像点: **条件成立, 必须外部验证**(见结论 1)。
%     * 径向尺度 f: 本方法**给不出**(见结论 2)。r 的无假设信息只有 R_rim 一个数。
%     * 畸变律 r(z): 不给(需要官方像素视角表或真实星场)。
%     * 方位角 rot: 任何"无星"方法都给不出 —— 本批数据的点源不随天空转。
%
%   见 mooncal_work\pipe_check\倾斜与天顶前提_检验.txt 与
%      mooncal_work\pipe_check\参数不确定度与退化_检验.txt (实测数字)。

%% ---------- argument parsing ----------
opt = struct('DarkFile', '', 'NFrames', 30, 'FovHalfDeg', NaN, ...
             'HalfFovRange', [70 110], 'ExportScale', false, 'ZenithRef', [], ...
             'ScalePxRad', NaN, ...
             'SearchHalf', 100, 'NSig', 2.5, 'SaveMat', '', 'SavePreview', '');
if mod(numel(varargin), 2) ~= 0
    error('zenith_from_rim:args', 'Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    nm = varargin{k};
    if ~isfield(opt, nm)
        error('zenith_from_rim:args', 'Unknown option "%s".', nm);
    end
    opt.(nm) = varargin{k+1};
end
opt.NFrames = max(1, min(60, round(opt.NFrames)));
opt.ExportScale = logical(opt.ExportScale);
fovKnown = isfinite(opt.FovHalfDeg) && opt.FovHalfDeg > 0;

if ischar(nightDir) || isstring(nightDir)
    dirs = {char(nightDir)};
else
    dirs = nightDir;
end

%% ---------- collect frames ----------
files = {};
for d = 1:numel(dirs)
    dd = dir(fullfile(dirs{d}, '*.PNG'));
    if isempty(dd), dd = dir(fullfile(dirs{d}, '*.png')); end
    for j = 1:numel(dd)
        files{end+1} = fullfile(dd(j).folder, dd(j).name); %#ok<AGROW>
    end
end
files = sort(files);
nAll = numel(files);
if nAll == 0
    error('zenith_from_rim:data', 'No PNG frames found under the given folder(s).');
end
step = max(1, floor(nAll / opt.NFrames));
sel  = files(1:step:end);
sel  = sel(1:min(numel(sel), opt.NFrames));
fprintf('>>> found %d frames, using %d\n', nAll, numel(sel));

%% ---------- dark frame ----------
darkImg = [];
if ~isempty(opt.DarkFile)
    if ~exist(opt.DarkFile, 'file')
        error('zenith_from_rim:dark', 'Dark file not found: %s', opt.DarkFile);
    end
    darkImg = q_read_gray(opt.DarkFile);   % 按位深统一到 16 bit 量级
    fprintf('>>> dark loaded: %s\n', opt.DarkFile);
end

%% ---------- median stack ----------
A = [];
for k = 1:numel(sel)
    a = q_read_gray(sel{k});            % 按位深统一到 16 bit 量级
    if isempty(a)
        warning('zenith_from_rim:decode', '解码失败，跳过: %s', sel{k});
        continue
    end
    if ~isempty(darkImg) && isequal(size(a), size(darkImg))
        a = a - darkImg;
    end
    if isempty(A)
        A = zeros([size(a), numel(sel)], 'single');
    end
    A(:, :, k) = single(a); %#ok<AGROW>
end
if isempty(A)
    error('zenith_from_rim:allBad', '所有帧都解码失败，无法做中值堆叠。');
end
med = double(median(A, 3, 'omitnan'));
[H, W] = size(med);
fprintf('>>> median stack done: %d x %d\n', W, H);

%% ---------- stage 1: half-maximum crossing, wide band -> robust coarse circle ----------
c0   = [W/2, H/2];
cNow = c0;
rg   = guess_radius(med, c0);
fprintf('>>> field radius (initial guess) %.1f px\n', rg);
cxB = NaN; cyB = NaN; rB = NaN; rmsB = NaN; nB = 0;
vIn = NaN; vOut = NaN;
for it = 1:3
    [P2, vIn, vOut] = edge_points_halfmax(med, cNow, rg, opt.SearchHalf);
    if size(P2, 1) < 60
        error('zenith_from_rim:edge', 'Too few half-max edge samples (%d).', size(P2, 1));
    end
    [cxB, cyB, rB, rmsB, nB] = fit_circle_robust(P2, opt.NSig);
    cNow = [cxB, cyB];
    rg   = rB;
end

%% ---------- stage 2: gradient peak, narrow band -> sub-pixel refinement ----------
% A narrow band is essential: a wide band lets the argmax jump to unrelated
% bright structure inside the field (bright fixed-pattern dots, cloud edges).
cNow = [cxB, cyB];
rg   = rB;
narrowHalf = 12;
cxA = NaN; cyA = NaN; rA = NaN; rmsA = NaN; nA = 0;
for it = 1:3
    P = edge_points_gradient(med, cNow, rg, narrowHalf);
    if size(P, 1) < 60
        error('zenith_from_rim:edge', 'Too few edge samples (%d).', size(P, 1));
    end
    [cxA, cyA, rA, rmsA, nA] = fit_circle_robust(P, opt.NSig);
    cNow = [cxA, cyA];
    rg   = rA;
end

fprintf('\n================= RESULT =================\n');
fprintf('A radial-gradient max : centre (%9.4f, %9.4f)  R = %8.4f  rms = %.4f  n = %d\n', ...
        cxA, cyA, rA, rmsA, nA);
if isfinite(cxB)
    fprintf('B half-max crossing   : centre (%9.4f, %9.4f)  R = %8.4f  rms = %.4f  n = %d\n', ...
            cxB, cyB, rB, rmsB, nB);
    fprintf('A - B                 : dx = %+.3f  dy = %+.3f  dR = %+.3f px\n', ...
            cxA-cxB, cyA-cyB, rA-rB);
else
    fprintf('B: not enough edge points\n');
end
fprintf('sky level (inside) %.1f | blind level (outside) %.1f\n', vIn, vOut);
fprintf('\n*** the centre above is the BORESIGHT (optical axis x detector), ***\n');
fprintf('*** NOT automatically the zenith point.  See the header.        ***\n');

%% ---------- focal: only if explicitly asked ----------
focal = NaN;
if fovKnown && opt.ExportScale
    focal = rA / deg2rad(opt.FovHalfDeg);
    fprintf('\n[EXPORTED] assumed field-edge zenith angle = %.2f deg -> focal = R/theta = %.4f px/rad\n', ...
            opt.FovHalfDeg, focal);
    fprintf('           *** THIS IS NOT A MEASUREMENT ***  the edge angle is an assumption.\n');
elseif fovKnown
    f0 = rA / deg2rad(opt.FovHalfDeg);
    fprintf('\n[NOT EXPORTED] you gave FovHalfDeg = %.2f deg -> R/theta = %.4f px/rad,\n', ...
            opt.FovHalfDeg, f0);
    fprintf('           but ExportScale is false, so focal stays NaN and rz_poly is NaN.\n');
    fprintf('           Radial scale must come from the moon/star calibration.\n');
else
    fprintf('\n[NO SCALE] FovHalfDeg not given (NaN) -> focal = NaN.  R_rim = %.4f px is the\n', rA);
    fprintf('           only assumption-free radial number this method produces.\n');
end

% 半径 -> f 的取值范围 (theta 未知时的系统区间)
lo = opt.HalfFovRange(1); hi = opt.HalfFovRange(2);
fprintf('\nf = R/theta for theta_rim in [%.0f, %.0f] deg  =>  f in [%.2f, %.2f] px/rad  (%.2fx spread)\n', ...
        lo, hi, rA/deg2rad(hi), rA/deg2rad(lo), (rA/deg2rad(lo))/(rA/deg2rad(hi)));
for th = (10:10:round(hi))
    if th >= lo
        fprintf('      theta_rim = %3d deg -> f = %8.2f px/rad\n', th, rA/deg2rad(th));
    end
end
fprintf('   => a single radius cannot pin the scale better than this factor. Do not feed it\n');
fprintf('      into the calibration chain as a measurement (see v6 [finding B]).\n');

%% ---------- optional: is the centre == the zenith point? ----------
zenCheck = [];
if ~isempty(opt.ZenithRef)
    xr = opt.ZenithRef(1); yr = opt.ZenithRef(2);
    d  = hypot(cxA - xr, cyA - yr);
    rough = false;
    if isfinite(opt.ScalePxRad) && opt.ScalePxRad > 0
        fForAng = opt.ScalePxRad;          % 月亮/星定标值 —— 正确用法
    elseif isfinite(focal) && focal > 0
        fForAng = focal;                   % 用户显式假设的 f
        rough = true;
    else
        fForAng = rA / deg2rad(90);        % 兜底: 很粗, 仅为了给个量级
        rough = true;
    end
    ddeg = rad2deg(d / fForAng);
    if d < 5
        verdict = 'HOLDS  (centre == zenith within 5 px)';
    elseif d < 10
        verdict = 'WEAK    (5~10 px: quote the offset, do not call it the zenith)';
    else
        verdict = 'FAILS   (>10 px: the rim centre is NOT the zenith point -- usable as an initial value only)';
    end
    fprintf('\n================= ZENITH PREMISE CHECK =================\n');
    fprintf('  rim centre (boresight) : (%9.3f, %9.3f)\n', cxA, cyA);
    fprintf('  reference zenith point : (%9.3f, %9.3f)\n', xr, yr);
    fprintf('  offset                 : %.3f px  =  %.3f deg  (f = %.2f px/rad%s)\n', ...
            d, ddeg, fForAng, ternary_str(rough, ', ROUGH', ', from ScalePxRad'));
    fprintf('  verdict                : %s\n', verdict);
    fprintf('  NOTE: the offset IS the boresight-vs-zenith tilt (f*t). It is measured here,\n');
    fprintf('        not assumed. If it grows, the camera was re-levelled / re-mounted, or the\n');
    fprintf('        boresight really is off zenith -> re-do this check on a fresh night.\n');
    zenCheck = struct('ref', [xr yr], 'd_px', d, 'd_deg', ddeg, 'verdict', verdict);
end

%% ---------- assemble cal struct ----------
cal = struct();
cal.zenith            = [cxA, cyA];    % 1-based; == boresight (see header)
cal.zenithIsBoresight = true;
cal.zenithIsZenith    = false;         % 只能靠外部验证; 别默认 true
cal.boresightVerified = ~isempty(zenCheck) && zenCheck.d_px < 10;
cal.radius            = rA;
cal.focal             = focal;         % NaN 除非显式假设 + ExportScale
cal.focalIsMeasurement = false;
if fovKnown && opt.ExportScale
    cal.focalSource = sprintf('ASSUMED field-edge zenith angle = %.2f deg (NOT a measurement)', ...
                              opt.FovHalfDeg);
else
    cal.focalSource = 'none: field-edge zenith angle unknown; use the moon/star calibration';
end
cal.k1 = 0;  cal.k2 = 0;
cal.rot = NaN;                          % NOT determined by this method
if isfinite(focal)
    cal.rz_poly = [0, 0, 0, 0, 0, 0, focal*180/pi, 0];
else
    cal.rz_poly = NaN(1, 8);            % 不写尺度: 让消费方显式失败/改用真标定
end
cal.provenance = sprintf(['field-edge circle fit; FovHalfDeg=%s; ExportScale=%d; dark=%s; ' ...
                          'centre=boresight; scale not measured'], ...
                         num2str(opt.FovHalfDeg), opt.ExportScale, opt.DarkFile);
cal.rotKnown = false;

zen = struct();
zen.center        = [cxA, cyA];
zen.center0       = [cxA-1, cyA-1];
zen.radius        = rA;
zen.rms           = rmsA;
zen.nPoints       = nA;
zen.focal         = focal;
zen.focalFromAssumedFov = fovKnown;
zen.halfFovAssumed      = opt.FovHalfDeg;
zen.centerHalf    = [cxB, cyB];
zen.dRadius       = rB - rA;
zen.cal           = cal;
zen.zenithCheck   = zenCheck;
zen.rotKnown      = false;

fprintf('\n>>> NOTE: rot is NOT determined by this method (NaN).\n');
fprintf('>>> centre(=boresight) and R_rim are usable; the radial scale and the azimuth\n');
fprintf('    angle must come from the moon/star calibration.\n');

%% ---------- save ----------
if ~isempty(opt.SaveMat)
    save(opt.SaveMat, 'cal', 'zen');
    fprintf('>>> written %s\n', opt.SaveMat);
end
if ~isempty(opt.SavePreview)
    make_preview(med, [cxA, cyA], rA, [cxB, cyB], rB, opt.SavePreview);
    fprintf('>>> preview %s\n', opt.SavePreview);
end
end

% ======================================================================
function s = ternary_str(cond, a, b)
if cond, s = a; else, s = b; end
end

% ======================================================================
% to_double 已删除：读帧一律走 q_read_gray（它同时负责位深对齐）。

% ======================================================================
function rg = guess_radius(img, c0)
[H, W] = size(img);
[XX, YY] = meshgrid(1:W, 1:H);
rr = hypot(XX-c0(1), YY-c0(2));
rmax = min(max(rr(:)), 0.55*min(W, H));
bins = 0:6:rmax;
prof = nan(size(bins));
for i = 1:numel(bins)-1
    m = rr >= bins(i) & rr < bins(i+1);
    if nnz(m) > 50
        prof(i) = median(img(m));
    end
end
v = prof(isfinite(prof));
if numel(v) < 5
    rg = 0.5*min(W, H);
    return;
end
inL  = v(1:max(3, round(0.35*numel(v))));
outL = v(max(1, round(0.88*numel(v))):end);
lv = median(inL);
lo = median(outL);
if lv - lo < 30
    rg = 0.47*min(W, H);
    return;
end
thr = lo + 0.35*(lv - lo);
idx = find(prof < thr, 1, 'first');
if isempty(idx)
    rg = 0.47*min(W, H);
else
    rg = bins(max(idx-1, 1));
end
rg = min(max(rg, 100), 0.55*min(W, H));
end

% ======================================================================
function P = edge_points_gradient(img, c0, rg, searchHalf)
[H, W] = size(img);
nang = 1440;
th = linspace(0, 2*pi, nang+1); th(end) = [];
rs = (rg-searchHalf):0.5:(rg+searchHalf);
X = c0(1) + cos(th(:))*rs;
Y = c0(2) + sin(th(:))*rs;
X(X < 1 | X > W) = NaN;
Y(Y < 1 | Y > H) = NaN;
V = interp2(img, X, Y, 'linear', NaN);
% smooth along the ray before differentiating: a plain 1-sample difference on
% 0.5 px spacing is dominated by read noise and locks onto the wrong feature.
% NOTE: do NOT use conv2(...,'same') here -- it zero-pads the ends and fabricates
% a huge artificial gradient at both band edges.  Clamp-replicate instead.
K  = [1 4 6 4 1]/16;
Vs = zeros(size(V));
nK = numel(K);
for q = 1:nK
    off = q - ceil(nK/2);
    ic  = min(max((1:size(V,2)) + off, 1), size(V,2));
    Vs  = Vs + K(q)*V(:, ic);
end
dV = diff(Vs, 1, 2) / 0.5;
absd = abs(dV);
[~, jj] = max(absd, [], 2);
nrs = size(absd, 2);
ok = isfinite(jj) & jj > 1 & jj < nrs;
jsub = double(jj);
for k = find(ok).'
    y0 = absd(k, jj(k)-1);
    y1 = absd(k, jj(k));
    y2 = absd(k, jj(k)+1);
    den = y0 - 2*y1 + y2;
    if abs(den) > 1e-12
        dd = 0.5*(y0 - y2)/den;
        jsub(k) = jj(k) + max(min(dd, 0.5), -0.5);
    end
end
pw = zeros(nang, 1);
pw(ok) = absd(sub2ind(size(absd), find(ok), jj(ok)));
ok = ok & isfinite(jsub) & (pw > 0.30*prctile_local(pw(ok), 90));
rEdge = rs(1) + (jsub - 1)*0.5;
thOk = th(ok).';
rOk  = rEdge(ok);
P = [c0(1) + rOk.*cos(thOk), c0(2) + rOk.*sin(thOk)];
end

% ======================================================================
function [P, vIn, vOut] = edge_points_halfmax(img, c0, rg, searchHalf)
[H, W] = size(img);
nang = 1440;
th = linspace(0, 2*pi, nang+1); th(end) = [];
rs = (rg-searchHalf):0.25:(rg+searchHalf);
X = c0(1) + cos(th(:))*rs;
Y = c0(2) + sin(th(:))*rs;
X(X < 1 | X > W) = NaN;
Y(Y < 1 | Y > H) = NaN;
V = interp2(img, X, Y, 'linear', NaN);
P = zeros(0, 2);
inAll = zeros(0, 1);
outAll = zeros(0, 1);
nb = min(40, floor(size(V, 2)/4));
for k = 1:nang
    p = V(k, :);
    if nnz(isfinite(p)) < 2*nb, continue; end
    vi = median(p(1:nb), 'omitnan');
    vo = median(p(end-nb+1:end), 'omitnan');
    if ~isfinite(vi) || ~isfinite(vo) || (vi - vo) < 50, continue; end
    half = (vi + vo)/2;
    j = find(p < half, 1, 'first');
    if isempty(j) || j == 1 || ~isfinite(p(j-1)), continue; end
    v0 = p(j-1); v1 = p(j);
    if abs(v1 - v0) < 1e-9, continue; end
    t = (half - v0)/(v1 - v0);
    rr = rs(j-1) + t*0.25;
    P(end+1, :) = [c0(1) + rr*cos(th(k)), c0(2) + rr*sin(th(k))]; %#ok<AGROW>
    inAll(end+1, 1)  = vi;  %#ok<AGROW>
    outAll(end+1, 1) = vo;  %#ok<AGROW>
end
if isempty(inAll)
    vIn = NaN; vOut = NaN;
else
    vIn = median(inAll); vOut = median(outAll);
end
end

% ======================================================================
function [cx, cy, r, rms, n] = fit_circle_robust(P, nsig)
q = P;
cx = NaN; cy = NaN; r = NaN; %#ok<NASGU>
for k = 1:6
    [cx, cy, r] = fit_circle_algebraic(q);
    d = hypot(q(:, 1)-cx, q(:, 2)-cy) - r;
    m = median(d);
    s = 1.4826*median(abs(d - m));
    if s <= 1e-9, break; end
    keep = abs(d - m) < nsig*s;
    if nnz(keep) < 40, break; end
    q = q(keep, :);
end
[cx, cy, r] = fit_circle_algebraic(q);
d = hypot(q(:, 1)-cx, q(:, 2)-cy) - r;
rms = sqrt(mean(d.^2));
n = size(q, 1);
end

% ======================================================================
function [cx, cy, r] = fit_circle_algebraic(P)
x = P(:, 1); y = P(:, 2);
AA = [2*x, 2*y, ones(size(x))];
bb = x.^2 + y.^2;
s = AA \ bb;
cx = s(1); cy = s(2);
r = sqrt(s(3) + cx^2 + cy^2);
end

% ======================================================================
function v = prctile_local(x, p)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
i = max(1, min(numel(x), round(p/100*numel(x))));
v = x(i);
end

% ======================================================================
function make_preview(img, cA, rA, cB, rB, fname)
f = figure('Color', 'w', 'Position', [100 100 900 860], 'Visible', 'off');
ax = axes('Parent', f);
v = img(:);
lo = prctile_local(v, 5);
hi = prctile_local(v, 99.5);
imagesc(ax, img, [lo, hi]);
colormap(ax, gray);
axis(ax, 'image');
hold(ax, 'on');
th = linspace(0, 2*pi, 720);
plot(ax, cA(1)+rA*cos(th), cA(2)+rA*sin(th), '-', 'Color', [1 0 0], 'LineWidth', 1.2);
plot(ax, cA(1), cA(2), '+', 'Color', [1 1 0], 'MarkerSize', 16, 'LineWidth', 2);
if isfinite(cB(1))
    plot(ax, cB(1)+rB*cos(th), cB(2)+rB*sin(th), '-', 'Color', [0 0.7 1], 'LineWidth', 1.0);
end
title(ax, sprintf('field edge circle   centre (%.2f, %.2f)   R = %.2f px', cA(1), cA(2), rA));
set(ax, 'FontSize', 10);
print(f, fname, '-dpng', '-r100');
close(f);
end
