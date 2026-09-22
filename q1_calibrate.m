function out = q1_calibrate(varargin)
%Q1_CALIBRATE  ① 定标一条龙：月亮法定标 → 写管线格式 cal → 几何法独立复核
%
%   q1_calibrate                       % 全默认，跑全部观测夜
%   q1_calibrate('Every', 2)           % 每 2 帧取 1，快一倍
%   q1_calibrate('Nights', {'20260322','20260323','20260324','20260325','20260326','20260327'})
%   q1_calibrate('SkipRim', true)      % 跳过几何复核
%
% 什么时候跑
%   ★ 同一台相机、没动过安装 ⇒ **不用重跑**，直接用现成的 calibration_params.mat。
%   ★ 相机拆装/重新调平/换站 ⇒ 必须重跑。
%   判据：跑一次 q4_verify 看几何圆心与老值 (543.82, 474.35) 差多少，
%         差 <2 px 属前者，差几十 px 属后者。
%
% 做什么（三步，原理见各函数文件头）
%   ① moon_calib.m          月亮法定出 (x0,y0) / rot / r(z)=f(z+a z^3)，最小二乘
%   ② moon_calib_to_pipeline 0-based→1-based、rot 符号不换、径向多项式转"度"制 8 元
%                            （由 moon_calib 的 'pipeline' 选项自动完成，并自动备份旧文件）
%   ③ zenith_from_rim.m     几何法独立测 boresight 与视场半径（不需要星、不需要月亮），
%                            与月亮天顶像点的差值就是"光轴偏天顶"的硬上限
%
% 参数（名值对）
%   'Nights'    参与定标的夜，默认 {} = 自动挑月亮合适的夜
%   'Every'     抽帧步长，默认 1
%   'Stage'     moon_calib 的阶段，默认 'all'
%   'MinIllum'  最低照亮比，默认 0.15
%   'AltLo'/'AltHi'   月亮高度窗口(deg)，默认 14 / 78
%   'SkipRim'   跳过 ③，默认 false
%   'RimNights' 几何复核用哪些夜，默认 {'20260516'}
%   'Chatty'    把 moon_calib 的详细长表也打到控制台，默认 false（只进日志文件）
%   'Verbose'   默认 true

p = inputParser;
addParameter(p, 'Nights', {});
addParameter(p, 'Every', 1);
addParameter(p, 'Stage', 'all');
addParameter(p, 'MinIllum', 0.15);
addParameter(p, 'AltLo', 14.0);
addParameter(p, 'AltHi', 78.0);
addParameter(p, 'SkipRim', false);
addParameter(p, 'RimNights', {'20260516'});
addParameter(p, 'Chatty', false);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

cfg = q_cfg();
[pn, bn] = fileparts(cfg.calNative);
% 原生定标（calMat）的**文件名**由 cfg.calNative 决定：用 'outmat' 直接落盘到该完整路径。
% 其余产物（calib_report.txt / json / csv / qc 图）落在同名子目录里 —— 这样把 rawRoot
% 切到另一批数据（例如一期 2014 的 8 bit 数据）时，两批定标的报告不会互相覆盖。
% ⚠ 历史坑：moon_calib 原生 .mat 的文件名**写死**为 calibration_params.mat；
%   只给 'out' 的话，第二批数据会静默覆盖第一批的定标文件。
if strcmp(bn, 'calibration_params')
    nativeOut = pn;
else
    nativeOut = fullfile(pn, bn);
end
if ~exist(nativeOut, 'dir'), mkdir(nativeOut); end
if ~exist(pn, 'dir'), mkdir(pn); end
if ~exist(fileparts(cfg.calFile), 'dir'), mkdir(fileparts(cfg.calFile)); end

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  ① 定标 · 月亮法\n');
fprintf('  原始数据 %s\n', cfg.rawRoot);
fprintf('  管线定标 %s\n', cfg.calFile);
fprintf('  原生定标 %s\n', cfg.calNative);
fprintf('%s\n', repmat('=', 1, 78));

args = {'raw', cfg.rawRoot, 'out', nativeOut, 'outmat', cfg.calNative, ...
    'site', 'official', ...
    'pipeline', cfg.calFile, 'stage', o.Stage, 'every', o.Every, ...
    'min_illum', o.MinIllum, 'alt_lo', o.AltLo, 'alt_hi', o.AltHi, ...
    'nights', o.Nights, 'chatty', o.Chatty};
moon_calib(args{:});

out = struct();
out.nativeMat = cfg.calNative;
out.pipeMat = cfg.calFile;

% ---------- ③ 几何法复核：光轴真的朝天顶吗 ----------
if ~o.SkipRim
    fprintf('\n%s\n', repmat('=', 1, 78));
    fprintf('  ③ 几何复核 · 视场边界圆拟合（zenith_from_rim）\n');
    fprintf('%s\n', repmat('-', 1, 78));
    fprintf('  原理：r 只是"离轴角"的函数 ⇒ 等天顶角的像点落在**以光轴像点为圆心**的圆上；\n');
    fprintf('        视场边界就是这样一条等角线 ⇒ 拟合边界圆的**圆心 = 光轴（boresight）像点**。\n');
    fprintf('  能力边界：它给不出 rot，也给不出 f ——\n');
    fprintf('        · 圆是旋转对称的，绕圆心转任意角度都不改变它 ⇒ 定不了方位零点；\n');
    fprintf('        · R_rim 对应多少天顶角需要外部参照 ⇒ f = R_rim/θ 里的 θ 图像自证不了。\n');
    fprintf('  所以本节只当**独立旁证**：它与月亮天顶像点之差 = "光轴偏天顶"的硬上限。\n');

    S = load(cfg.calFile);
    zref = [double(S.cal.zenith(1)), double(S.cal.zenith(2))];   % 1-based
    nd = o.RimNights;
    if ischar(nd), nd = {nd}; end
    dirs = cellfun(@(n) fullfile(cfg.rawRoot, n), nd, 'UniformOutput', false);

    z = zenith_from_rim(dirs, 'NFrames', 12, 'ZenithRef', zref, ...
        'ScalePxRad', double(S.cal.rz_poly(end - 1)) * 180 / pi);
    dpx = hypot(z.center(1) - zref(1), z.center(2) - zref(2));
    fp = double(S.cal.rz_poly(end - 1)) * 180 / pi;
    fprintf('%s\n', repmat('-', 1, 78));
    fprintf('  几何法圆心（光轴像点） %10.3f, %10.3f    1-based\n', z.center(1), z.center(2));
    fprintf('  月亮法天顶像点         %10.3f, %10.3f\n', zref(1), zref(2));
    fprintf('  边界圆半径 R_rim       %10.3f px   （拟合 rms %.4f px）\n', z.radius, z.rms);
    fprintf('  两者相差               %10.2f px = %.3f°\n', dpx, atand(dpx / fp));
    fprintf('  判定：');
    if dpx < 10
        fprintf('一致（< 10 px）⇒ "天顶像点 = 光轴像点"的前提成立。\n');
    else
        fprintf('偏差偏大 ⇒ 检查是否重新安装/调平过相机，或换过站。\n');
    end
    fprintf('%s\n', repmat('=', 1, 78));
    out.rim = z;
end
end
