function out = start_calib(varargin)
%START_CALIB  ★启动脚本 ①  定标★
%
%   start_calib                        % 全默认：自动挑月亮合适的夜，跑完整定标
%   start_calib('Every', 2)            % 每 2 帧取 1，快一倍
%   start_calib('Nights', {'20260326','20260327'})
%   start_calib('SkipRim', true)       % 跳过几何复核
%   start_calib('Verify', false)       % 定标后不跑自检
%   start_calib('DryRun', true)        % 只看环境与现有定标，不做任何计算与写入
%
% ■ 什么时候需要跑
%   同一台相机、没动过安装   ⇒ **不用重跑**，直接用现成的 calibration_params.mat
%   相机拆装 / 重新调平 / 换站 ⇒ **必须重跑**
%   判别办法：跑一次 q4_verify，看"几何圆心"与"月亮天顶像点"差多少 ——
%             差 < 2 px 属前者（安装没变），差几十 px 属后者（必须重跑）。
%
% ■ 本脚本做三件事（原理见各函数文件头）
%   ① moon_calib.m            月亮法定出 (x0,y0) / rot / r(z) = f (z + a z^3)，最小二乘
%   ② moon_calib_to_pipeline  0-based→1-based、rot 符号照抄、径向多项式转"度"制 8 元，
%                             自动备份旧文件（由 moon_calib 的 'pipeline' 选项完成）
%   ③ zenith_from_rim.m       几何法独立测光轴像点与视场半径。它给不出 rot 和 f，
%                             但与月亮天顶像点之差 = "光轴偏天顶"的硬上限
%                             —— 这是**测出来的**，不是假设的。
%
% ■ 产出
%   <cfg.calFile>                            管线用定标（cal 结构：zenith/rot/rz_poly/...）
%   <cfg.calNative>                          moon_calib 原生定标（calMat，0-based）
%   <cfg.logRoot>\start_calib_<时间>.txt     本次全过程日志
%
% ■ 换原始图像文件夹
%   只改 q_cfg.m 第 21 行 cfg.rawRoot。本脚本从 q_cfg 取全部路径，别处不用动。
%
% ■ 参数（名值对）
%   'Nights'  参与定标的夜，默认 {} = 自动挑月亮合适的夜
%   'Every'   抽帧步长，默认 1
%   'SkipRim' 跳过 ③ 几何复核，默认 false
%   'Verify'  定标后跑 q4_verify 自检，默认 true
%   'DryRun'  只回显环境与现有定标，不计算、不写文件，默认 false
%   'Chatty'  把月亮法的详细长表也打到控制台，默认 false（长表只进日志文件）
%   'Diary'   全过程写日志文件，默认 true
%   'Verbose' 默认 true
%
% ■ 输出详略
%   默认**精简**：控制台只打结论与判读，末尾给一张「定标结果卡」；逐夜/逐帧长表全量写进
%   <cfg.logRoot>\start_calib_<时间>.txt，事后可复算。加 'Chatty',true 可现场看全量。
%
% 输出 out: nativeMat / pipeMat / cal / rim / verify / logFile / elapsedSec

p = inputParser;
addParameter(p, 'Nights', {});
addParameter(p, 'Every', 1);
addParameter(p, 'SkipRim', false);
addParameter(p, 'Verify', true);
addParameter(p, 'DryRun', false);
addParameter(p, 'Chatty', false);
addParameter(p, 'Diary', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
o = p.Results;

t0 = tic;

% ==========================================================================
% 0  定位工程根（脚本在哪儿，工程根就在哪儿 —— 不写死盘符）
% ==========================================================================
here = fileparts(mfilename('fullpath'));
if ~isempty(here)
    cd(here);
    addpath(here);
end
if exist('q_cfg.m', 'file') ~= 2
    error('start_calib:noCfg', '找不到 q_cfg.m —— 本脚本必须与它放在同一目录。');
end

% ==========================================================================
% 0.1  日志（diary 落盘；函数退出/报错时自动关闭）
% ==========================================================================
cfg = q_cfg();
logFile = '';
if o.Diary
    if exist(cfg.logRoot, 'dir') ~= 7, mkdir(cfg.logRoot); end
    logFile = fullfile(cfg.logRoot, ...
        ['start_calib_' datestr(now, 'yyyymmdd_HHMMSS') '.txt']);
    diary(logFile); diary on;
    % 这个变量必须存活到函数结束（含中途报错），diary 才会被关掉
    cleanupDiary = onCleanup(@() diary('off'));
end

% ==========================================================================
fprintf('\n%s\n', repmat('#', 1, 78));
fprintf('#  启动脚本 ①  定标   start_calib\n');
fprintf('#  时间 %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('#  工程根 %s\n', here);
fprintf('%s\n', repmat('#', 1, 78));

sec_hdr = @(s) fprintf('\n%s\n  %s\n%s\n', repmat('=', 1, 78), s, repmat('=', 1, 78));

sec_hdr('0  路径体检');
q_cfg('show');
q_cfg('check');

% ==========================================================================
sec_hdr('0.1  是否需要重定标');
fprintf('  同相机、没动过安装        ⇒ 不必重跑，直接用现成的 .mat\n');
fprintf('  相机拆装/重新调平/换站    ⇒ 必须重跑（本脚本就干这个）\n');
fprintf('  判别：q4_verify 里"几何圆心"与"月亮天顶像点"之差 < 2 px 说明安装没变。\n');

sec_hdr('0.2  现有定标（本次将覆盖）');
if exist(cfg.calFile, 'file') == 2
    S = load(cfg.calFile);
    if isfield(S, 'cal')
        show_cal(S.cal);
    else
        fprintf('  文件里没有 cal 变量：%s\n', cfg.calFile);
    end
else
    fprintf('  尚无管线定标文件（首次定标）：%s\n', cfg.calFile);
end

% ==========================================================================
% 1  月亮法定标（内部：moon_calib → to_pipeline → 几何复核）
% ==========================================================================
r1 = [];
if o.DryRun
    fprintf('\n%s\n', repmat('=', 1, 78));
    fprintf('  ① 定标（DryRun：已跳过）\n');
    fprintf('  只做环境体检与现有定标回显，不计算、不写任何文件。\n');
    fprintf('%s\n', repmat('=', 1, 78));
else
    r1 = q1_calibrate('Nights', o.Nights, 'Every', o.Every, ...
        'SkipRim', o.SkipRim, 'Chatty', o.Chatty, 'Verbose', o.Verbose);
end

% ==========================================================================
% 2  读回定标，打印摘要
% ==========================================================================
out = struct();
out.nativeMat = cfg.calNative;
out.pipeMat = cfg.calFile;
out.cal = [];
out.rim = [];
out.verify = [];
if ~isempty(r1) && isfield(r1, 'rim'), out.rim = r1.rim; end

sec_hdr('2  定标结果');
if exist(cfg.calFile, 'file') == 2
    S = load(cfg.calFile);
    if isfield(S, 'cal')
        out.cal = S.cal;
        show_cal(S.cal);
    elseif o.DryRun
        fprintf('  [警告] 文件里没有 cal 变量：%s\n', cfg.calFile);
    else
        error('start_calib:badOut', '定标文件缺少 cal 变量：%s', cfg.calFile);
    end
elseif o.DryRun
    fprintf('  [警告] 尚无管线定标文件：%s（首次使用请去掉 DryRun 跑一次真定标）\n', cfg.calFile);
else
    error('start_calib:noOut', '定标跑完却没有产出管线定标：%s', cfg.calFile);
end

% 交叉核对：原生定标是否也写出来了
if exist(cfg.calNative, 'file') == 2
    fprintf('\n  原生定标已写出：%s\n', cfg.calNative);
else
    fprintf('\n  [警告] 原生定标未见：%s（自检 ⑤ 会用到）\n', cfg.calNative);
end

% ==========================================================================
% 3  定标后自检
% ==========================================================================
if o.Verify && ~o.DryRun
    sec_hdr('3  定标后自检 q4_verify');
    try
        out.verify = q4_verify('SkipSlow', true, 'Verbose', o.Verbose);
    catch ME
        fprintf('  [自检失败] %s\n', ME.message);
    end
else
    fprintf('\n（已跳过定标后自检）\n');
end

% ==========================================================================
out.logFile = logFile;
out.elapsedSec = toc(t0);

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('定标流程结束，用时 %.1f s\n', out.elapsedSec);
fprintf('  管线定标 : %s\n', cfg.calFile);
fprintf('  原生定标 : %s\n', cfg.calNative);
if ~isempty(logFile)
    fprintf('  日志     : %s\n', logFile);
end
fprintf('下一步   : start_process(''20250920'')   做图像处理 + 地理投影\n');
fprintf('%s\n', repmat('=', 1, 78));
end

% ==========================================================================
function show_cal(cal)
%SHOW_CAL 定标摘要卡（管线用 cal 结构）。
%
% 为什么不再逐字段刷屏：原实现把 cal 的 ~20 个字段全打出来，其中大半是约定说明
% （rot_convention / rz_poly_model / source …），每天真正要看的只有 4 项 ——
% 天顶像点、rot、径向标度、站点。这里只留这些 + 两条易错提示，全字段仍在 .mat 里。
W = 74;
fprintf('  %s\n', repmat('-', 1, W));
if isfield(cal, 'zenith')
    fprintf('  天顶像点（画布原点）  (%9.3f, %9.3f)   [1-based]\n', cal.zenith(1), cal.zenith(2));
    fprintf('                        0-based 时相当于 (%9.3f, %9.3f)\n', ...
            cal.zenith(1) - 1, cal.zenith(2) - 1);
end
if isfield(cal, 'rot')
    fprintf('  方位校正角 rot        %+.4f°   （投影由逆映射一次做掉，不用单独旋转）\n', cal.rot);
end
if isfield(cal, 'rz_poly')
    q = double(cal.rz_poly(:))';
    if numel(q) >= 2
        fpx = q(end - 1) * 180 / pi;
        fprintf('  径向标度 r(z)         r(90°) = %.2f px   ⇒ f ≈ %.2f px/rad（%.4f px/deg）\n', ...
                polyval(q, 90), fpx, fpx * pi / 180);
    end
end
if isfield(cal, 'radius')
    fprintf('  视场半径 radius       %.2f px\n', cal.radius);
end
if isfield(cal, 'z_rim_deg') && ~isempty(cal.z_rim_deg) && isfinite(cal.z_rim_deg)
    if cal.z_rim_deg < 95
        fprintf('  盘边天顶角 z_rim      %.2f°   （**不是 90°** ⇒ 勿用 f = R_rim/(π/2) 反推）\n', ...
                cal.z_rim_deg);
    else
        fprintf('  盘边天顶角 z_rim      %.2f°   （> 95° 无物理含义：radius 只是画幅裁边，非光学盘边）\n', ...
                cal.z_rim_deg);
    end
end
if isfield(cal, 'obs_lat') && isfield(cal, 'obs_lon')
    fprintf('  站点                  %.4f°N, %.4f°E\n', cal.obs_lat, cal.obs_lon);
end
if isfield(cal, 'rms_px')
    fprintf('  拟合质量              rms %.3f px', cal.rms_px);
    if isfield(cal, 'n_frames') && isfield(cal, 'n_nights')
        fprintf('（%d 帧 / %d 夜）', cal.n_frames, cal.n_nights);
    end
    fprintf('\n');
end
if isfield(cal, 'focal') && isfield(cal, 'a')
    fprintf('  径向模型              r = %.6g·(z + %+.6g·z³)，z 用弧度\n', cal.focal, cal.a);
end
fprintf('  %s\n', repmat('-', 1, W));
fprintf('  全字段：load 定标 .mat 后 disp(cal)。本卡不逐字段罗列，字段增删无需改这里。\n');
end
