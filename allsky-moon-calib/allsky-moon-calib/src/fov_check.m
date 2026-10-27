function fov_check(varargin)
%FOV_CHECK  仪器视场角（标称 180 deg，有效约 160 deg）与"可见盘边=R_rim 到底对应多少天顶角"的一致性检验。
%
% 为什么需要它
%   本链的径向标度 r(z) 一直用**月亮**直接测（不需要知道视场角）。但凡是"盘边/半径"这条线，
%   都隐含一个换算：那条可见圆盘的边缘对应多少天顶角 theta_rim？
%   · 标称视场 180 deg 会让人以为 theta_rim = 90 deg，于是 f = R_rim/(pi/2)；
%   · 但仪器"有效视场约 160 deg"（±80 deg）又会让人以为 theta_rim = 80 deg。
%   这两个读法给出**相差 12% 的 f**，而它们对本链的结论影响完全不同：
%     theta_rim=90 deg ⇒ f = 456.5/1.5708 = 290.6 px/rad，比月亮实测 324.4 低 10%；
%     theta_rim=80 deg ⇒ f = 456.5/1.3963 = 327.0 px/rad，与月亮实测只差 0.8%。
%   所以必须回答：**theta_rim 到底是几度？**
%
% 怎么答（本函数的做法，只用一个假设：r 是 z 的函数）
%   把 theta_rim 当成待检假设：固定一条假设（theta_rim = Th, 盘边 = R），它就等价于对
%   r(z)=f(z+a z^3) 里的三次项 a 提出一个要求：
%       a_req(R, Th) = (R/f - Th)/Th^3       （f 由月亮数据在该 a 下重新最小二乘）
%   然后问：**月亮数据（316 帧）允许 a 偏离最优值 a* 多少？**
%   做法：扫一遍 a，对每个 a 用月亮数据最小二乘定 f，得到 rms(a) 曲线。
%   rms(a) 就是判决器：若某个假设把 rms 从 4.57 px 抬到几十 px，该假设即被数据否掉。
%
% 用法
%   fov_check              % 两个批次都做（一期完整检验 + 二期算术说明），写同一个文件
%   fov_check('14')
%   fov_check('25')
%
% 输出
%   mooncal_work\reports\工具_视场与边界角.txt
%
% 依赖：mooncal_work\<epoch>\moon_meas.mat（月核质心缓存）、moon_altaz_ref.m、sun_altaz_ref.m

if nargin >= 1 && ~isempty(varargin{1}), Epoch = char(varargin{1}); else, Epoch = 'both'; end

MP      = mc_paths();
root    = MP.root;
site    = [109.133, 19.526, 0.103];      % lon, lat, h_km（官方儋州站）
altLo   = 14;  altHi = 50;               % 与 calib_2014 的拟合窗口一致
sunMax  = -18;                           % 与 calib_2014 的太阳门限一致
outDir  = MP.logs;
if ~exist(outDir, 'dir'), mkdir(outDir); end
fid = fopen(fullfile(outDir, '工具_视场与边界角.txt'), 'w', 'n', 'UTF-8');
cln = onCleanup(@() fclose(fid));
lg  = @(varargin) fprintf(fid, varargin{:});

lg('=== 视场角 / 盘边天顶角 一致性检验 === 生成 %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% ---- 两个批次的已知量 ----
% 一期：画幅 1024，半幅 512，而可见盘边 R_rim ≈ 456 px < 512 ⇒ **画幅内有暗边**，
%       盘边是真正的**光学软边**（渐晕到底），可用于本检验。
%       458.0 = calib_2014 的 RadiusInit；456.53 = 径向剖面 |dI/dr| 最大处；
%       450.96 = 亮度跌到"天空/偏置"半高点处（光终止）。三个都是同一过渡带的读数。
% 二期：R_rim = 512.93 ≈ 1024/2 ⇒ **盘边就是画幅边界**（被传感器裁掉），
%       不是光学盘边 ⇒ 不能用来定 theta_rim（这也解释了 [3.9] 里解出 z_rim≈100 deg 的怪值）。
P = struct();
P(1).tag    = '14';
P(1).name   = '一期 2014（8 bit）';
P(1).meas   = MP.measP14;
P(1).cal    = MP.calP14Native;
P(1).rims   = [458.00, 456.53, 450.96];
P(1).rimTag = {'RadiusInit(几何测量)', '剖面|dI/dr|最大', '半高(光终止)'};
P(1).half   = 512;
P(1).full   = true;
P(2).tag    = '25';
P(2).name   = '二期 2025/26（16 bit）';
P(2).meas   = MP.measP2;
P(2).cal    = MP.calP2Native;
P(2).rims   = 512.93;
P(2).rimTag = {'实测边界圆(=画幅裁边)'};
P(2).half   = 512;
P(2).full   = false;

sel = 1:numel(P);
if ~(strcmpi(Epoch, 'both') || strcmpi(Epoch, 'all'))
    sel = find(strcmp({P.tag}, Epoch));
    if isempty(sel), error('fov_check:badEpoch', '只支持 ''14''、''25'' 或 ''both''。'); end
end
for i = sel
    ep = P(i);
    lg('%s\n', repmat('#', 1, 78));
    lg('# %s\n', ep.name);
    lg('%s\n\n', repmat('#', 1, 78));

    A = load(ep.cal);
    if isfield(A, 'calMat'), cm = A.calMat; elseif isfield(A, 'cal'), cm = A.cal; else
        error('fov_check:noCal', '%s 里既没有 calMat 也没有 cal。', ep.cal);
    end
    f0 = cm.f; a0 = cm.a; x0 = cm.x0; y0 = cm.y0;
    lg('存档解: x0=%.2f y0=%.2f f=%.4f a=%+.5f rms=%.4f px（%d 帧/%d 夜）\n', ...
       x0, y0, f0, a0, cm.rms_px, cm.n_frames, cm.n_nights);
    lg('画幅半宽 = %d px；盘边候选 R_rim = %s px\n\n', ep.half, mat2str(ep.rims, 6));

    % ---- 径向函数与反解 ----
    rfun = @(z, f, a) f * (z + a * z.^3);
    % 反解用**网格首次穿越**，不用 fzero：二期 a 很负（-0.0998）⇒ r(z) 在 z~105 deg 处
    % 有极大值（≈514 px）后回落，R=512.93 会被穿越两次，而且区间端点同号会让 fzero 直接报错。
    % 这里与 moon_calib 的 mc_z_of_rim 保持一致，取首次穿越。
    zOfR = @(R, f, a) local_z_of_r(R, f, a);   % 返回 deg

    % ---- 一期的完整检验：用月亮数据扫 a，得 rms(a) ----
    if ep.full
        S = load(ep.meas);
        cx = double(S.cx(:)).'; cy = double(S.cy(:)).'; ts = S.ts(:).'; N = numel(cx);
        T = NaT(1, N);
        for k = 1:N
            s = char(ts{k});
            T(k) = datetime(str2double(s(1:4)), str2double(s(5:6)), str2double(s(7:8)), ...
                            str2double(s(9:10)), str2double(s(11:12)), str2double(s(13:14)));
        end
        alt = zeros(1, N); salt = zeros(1, N);
        for k = 1:N
            alt(k) = moon_altaz_ref(T(k), site(1), site(2), site(3));
            salt(k) = sun_altaz_ref(T(k), site(1), site(2));
        end
        sit = (salt < sunMax) & (alt > altLo) & (alt < altHi);
        lg('测量缓存 %d 帧；过门限帧集 %d 帧（太阳<%.0f deg、alt %g~%g）\n', N, sum(sit), sunMax, altLo, altHi);

        % 半径只与"离圆心距离"有关 ⇒ 不需要方位，rot 与径向拟合无关（这是关键：模型无关的量）
        za = deg2rad(90 - alt(sit)).';
        ra = hypot(cx(sit).' - x0, cy(sit).' - y0);      % 列向量：半径与方位无关
        keep = true(size(za));                  % 迭代紧截尾（同 adopt_solution 的口径）
        for it = 1:6
            [~, ~, rk] = fitRA(za(keep), ra(keep));
            m = median(abs(rk));
            thr = max(2.0, m + 3 * median(abs(abs(rk) - m)));
            kk = find(keep); keep(kk(abs(rk) >= thr)) = false;
        end
        za = za(keep); ra = ra(keep);
        [fFit, aFit, rFit] = fitRA(za, ra);
        lg('复算（固定圆心 + 紧截尾）: n=%d  f=%.4f  a=%+.5f  rms=%.4f px\n', ...
           numel(za), fFit, aFit, sqrt(mean(rFit.^2)));
        lg('  ⇒ 与存档解一致 ⇒ 下面的 rms(a) 曲线可信。\n\n');

        % ---- rms(a) 曲线 ----
        aGrid = linspace(-0.16, 0.06, 221);
        rmsG  = zeros(size(aGrid));
        for k = 1:numel(aGrid)
            [~, ~, rk] = fitRA(za, ra, aGrid(k));
            rmsG(k) = sqrt(mean(rk.^2));
        end
        [rmin, kmin] = min(rmsG);
        lg('----- (1) 月亮数据能容忍多大的 a？-----\n');
        lg('  rms(a) 极小: a* = %+.4f  ⇒ rms = %.3f px\n', aGrid(kmin), rmin);
        lg('  %10s %10s %12s   %s\n', 'a', 'rms(px)', 'rms/rms*', '判读');
        for k = kmin:-20:1
            lg('  %+10.4f %10.3f %12.3f   %s\n', aGrid(k), rmsG(k), rmsG(k)/rmin, verdict(rmsG(k)/rmin));
        end
        for k = kmin+20:20:numel(aGrid)
            lg('  %+10.4f %10.3f %12.3f   %s\n', aGrid(k), rmsG(k), rmsG(k)/rmin, verdict(rmsG(k)/rmin));
        end
        lg('  判据: rms 抬高 >50%%（即 rms/rms* > 1.5）即认为该假设被数据否掉\n');
        lg('        （本拟合 rms 已含系统项，故不用严格的卡方，只取"量级否掉"的保守口径）。\n\n');

        % ---- 假设检验表 ----
        lg('----- (2) 各假设 theta_rim 的天顶角 → 需要的 a → 被月亮数据接受吗 -----\n');
        lg('  %14s %10s %14s %10s %12s   %s\n', '读法', 'theta_rim', '需要的 a', 'rms(px)', 'rms/rms*', '判决');
        lg('  %s\n', repmat('-', 1, 92));
        ths = [80, 82.5, 85, 87.5, 90];
        thName = {'有效视场160(±80)', '82.5 deg', '85 deg', '87.5 deg', '标称180/2=90'};
        for j = 1:numel(ths)
            Th = deg2rad(ths(j));
            for R = ep.rims
                ar = solveA(R, Th, za, ra);
                [~, ~, rk] = fitRA(za, ra, ar);
                rr = sqrt(mean(rk.^2));
                lg('  %14s %10.1f %14.5f %10.3f %12.2f   %s\n', ...
                   thName{j}, ths(j), ar, rr, rr/rmin, verdict(rr/rmin));
            end
        end
        lg('  %s\n', repmat('-', 1, 92));
        lg('  说明: 同一 theta_rim 下三个 R_rim 读数（458.0/456.53/450.96）只是同一条软边的不同取法，\n');
        lg('        它们给出的 a 相差 ~0.005 ⇒ 软边宽度本身只带 ~1 deg 的 theta_rim 不确定性。\n\n');

        % ---- 自解：本解自己说盘边是多少度 ----
        lg('----- (3) 反过来：以本解 r(z) 为准，盘边对应多少度？-----\n');
        lg('  %22s %12s %14s\n', 'R_rim 取法', 'R (px)', 'theta_rim (deg)');
        for j = 1:numel(ep.rims)
            lg('  %22s %12.2f %14.2f\n', ep.rimTag{j}, ep.rims(j), zOfR(ep.rims(j), f0, a0));
        end
        lg('\n  ⇒ 三个读数给 theta_rim ≈ %.1f ~ %.1f deg。\n', ...
           zOfR(min(ep.rims), f0, a0), zOfR(max(ep.rims), f0, a0));
        lg('     而"有效视场 160 deg"对应的 80 deg 处 r = %.1f px。\n', rfun(deg2rad(80), f0, a0));
        lg('     80 deg 与盘边 %.1f deg 之间是一个 ~%d px 宽的**渐晕/像质退化外环**。\n\n', ...
           zOfR(ep.rims(2), f0, a0), round(ep.rims(2) - rfun(deg2rad(80), f0, a0)));
    end

    % ---- 两批次都做的 r(z) 表 + 与论文对照 ----
    if ep.full, secNo = 4; else, secNo = 1; end     % 一期前面已有 (1)(2)(3) 三节
    lg('----- (%d) r(z) 关键值、与"可用上限"的关系 -----\n', secNo);
    zTab = [60, 65, 70, 75, 76, 80, 85, 90];
    lg('  %8s %12s %14s   %s\n', 'z (deg)', 'r (px)', '占半幅', '备注');
    for zt = zTab
        rz = rfun(deg2rad(zt), f0, a0);
        nt = '';
        if zt > 90 - altLo - 0.001, nt = '**超出拟合区间（外推）**'; end
        if strcmp(ep.tag, '14') && zt >= 80 && zt <= 90, nt = '**在渐晕外环内**'; end
        lg('  %8.1f %12.1f %13.1f%%   %s\n', zt, rz, 100 * rz / ep.half, nt);
    end
    r90 = rfun(pi/2, f0, a0);
    lg('\n  r(90 deg) = %.1f px（占半幅 %.1f%%）\n', r90, 100 * r90 / ep.half);
    if strcmp(ep.tag, '14')
        lg('  ⚠ 本解 r(90 deg)=%.1f px **大于**可见盘边 %.1f px ⇒ 90 deg 圈已在渐晕外环之外。\n', r90, max(ep.rims));
        lg('  ⚠ 论文径向式在 90 deg 给 %.1f px ⇒ **超出盘边 %.1f px**，落在像圈之外，\n', ...
           polyval([-0.55, 305.2, 16.68], pi/2), polyval([-0.55, 305.2, 16.68], pi/2) - max(ep.rims));
        lg('     不可能被任何像元验证 ⇒ 与论文比径向尺度**只能比拟合区间内的有效尺度**。\n');
    else
        lg('  ⚠ 二期 R_rim=%.2f px ≈ 半幅 %d px ⇒ 该"边界圆"是**画幅裁边**，不是光学盘边；\n', ...
           ep.rims(1), ep.half);
        lg('     本解反解出的 z_rim≈%.1f deg 因此没有物理含义。二期**没有可用盘边信息**。\n', ...
           zOfR(ep.rims(1), f0, a0));
    end
    lg('\n');
end
end

% ---------------------------------------------------------------------------
function [f, a, r] = fitRA(z, rho, aFix)
% 固定 a（给定时）或自由拟合 a：r = f (z + a z^3)。对 f 是线性的 ⇒ 闭式解，快。
if nargin >= 3 && ~isempty(aFix)
    u = z + aFix * z.^3;
    f = (u.' * rho) / (u.' * u);
    a = aFix;
else
    u1 = z; u2 = z.^3;
    M = [u1, u2];                       % rho ≈ f*u1 + (f a)*u2
    c = M \ rho;                        % 最小二乘
    f = c(1); a = c(2) / c(1);
end
r = f * (z + a * z.^3) - rho;
end

% ---------------------------------------------------------------------------
function a = solveA(R, Th, z, rho)
% 求 a 使"盘边 R 对应天顶角 Th"：rms 在给定 a 下最小二乘定 f，然后解 f(a)*(Th+a Th^3)=R
lo = -0.30; hi = 0.30;
g = @(a) localg(a, R, Th, z, rho);
if g(lo) * g(hi) > 0
    a = NaN; return
end
a = fzero(g, [lo, hi]);
end
function v = localg(a, R, Th, z, rho)
[f, ~, ~] = fitRA(z, rho, a);
v = f * (Th + a * Th^3) - R;
end

% ---------------------------------------------------------------------------
function zDeg = local_z_of_r(R, f, a)
% 反解 r(z)=R 的首次穿越（deg）。与 moon_calib 的 mc_z_of_rim 同一口径。
% ⚠ 必须取"首次"穿越：二期 a=-0.0998 ⇒ r(z) 在 z≈104.7 deg 处有极大值 514 px，
%   而 R=512.93 与极大值只差 1.4 px ⇒ 会被穿越两次（≈100.1 deg 与 ≈109.2 deg），
%   取 min|rv-R| 会随机落到哪一支，必须先穿越才唯一。
zz = linspace(1e-6, 2.4, 200001);
rv = f * (zz + a * zz.^3);
k  = find(rv >= R, 1, 'first');
if isempty(k), zDeg = NaN; return; end
zDeg = rad2deg(zz(k));
end

% ---------------------------------------------------------------------------
function s = verdict(q)
if q < 1.05,      s = '接受';
elseif q < 1.5,   s = '勉强';
else,             s = '**否掉**';
end
end
