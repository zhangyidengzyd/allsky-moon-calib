function out = run_demo(mode, varargin)
%RUN_DEMO  ★整包唯一入口★
%
%   run_demo                     % = run_demo('info')：包概览 + 路径体检 + 定标结果卡
%   run_demo('check')            % 三项硬自检 + 论文数值逐条核对（**不需要原始数据**）
%   run_demo('project','20250920')   % 地理投影（需要该夜的原始帧，见 README「数据放置」）
%   run_demo('bubble','20250920')    % 泡团检测与漂移量测（需要该夜原始帧）
%   run_demo('calib','DryRun',true)  % 二期定标（先 DryRun 看环境，确认后再跑真的）
%   run_demo('calib','Nights',{'20260326','20260327'})
%   run_demo('phase1')               % 一期 2013/2014 定标（calib_2014）
%   run_demo('moon','20260326')      % 月夜专线：月面圆心归一化 + 月盘掩膜 + 探测极限
%
% ■ 这个包能做什么、不能做什么（先看这段再决定要不要下数据）
%   能：定标（月亮法）、几何投影、泡团检测与漂移量测、以及**全部论文数值的复算**。
%   不能：整链需要原始帧才能端到端复跑；原始帧共约 240 GB，**未随包提供**。
%         但 data/calib/ 里随包给了三个时期的采用解与逐帧月核量测缓存，
%         所以 §3 的定标结果与 §4 的投影/漂移数值**不需要原始数据**就能核对
%         （这正是 run_demo('check') 做的事）。
%
% ■ 想只用一分钟确认"这包是好的"
%   run_demo('check')
%
% 详见 README.md。

if nargin < 1 || isempty(mode), mode = 'info'; end
mode = lower(char(mode));

B = fileparts(mfilename('fullpath'));            % 包根
addpath(fullfile(B, 'src'));
addpath(fullfile(B, 'src', 'verify'));

out = struct('mode', mode, 'ok', true);

switch mode
% ---------------------------------------------------------------- info
case 'info'
    fprintf('\n%s\n', repmat('#', 1, 78));
    fprintf('#  全天空气辉成像仪 月亮法定标 + 气辉投影   代码包\n');
    fprintf('#  包根 %s\n', B);
    fprintf('%s\n', repmat('#', 1, 78));
    q_cfg('show');
    q_cfg('check');
    fprintf('\n--- 随包标定结果（摘要）---\n');
    show_card('二期 2026-03', mc_paths().calP2Native);
    show_card('一期 2014-09', mc_paths().calP14Native);
    show_card('一期 2013-10', mc_paths().calP13Native);
    fprintf('\n提示：run_demo(''check'') 跑三项硬自检 + 论文数值逐条核对。\n\n');

% ---------------------------------------------------------------- check
case 'check'
    a1 = verify_paths();
    a2 = verify_paper_numbers();
    out.paths = a1; out.numbers = a2;
    out.ok = (a1.bad == 0) && (a2.bad == 0) && (a2.miss == 0);

% ---------------------------------------------------------------- calib
case 'calib'
    % 二期：重新用月亮定标（会覆盖 data/calib/phase2_*.mat，先备份！）
    out.r = start_calib(varargin{:});

% ---------------------------------------------------------------- phase1
case 'phase1'
    % 一期 2013/2014：几何定标（月亮方位共识 + 固定圆心拟合尺度）
    out.r = calib_2014(varargin{:});

% ---------------------------------------------------------------- project
case 'project'
    if isempty(varargin)
        error('run_demo:needNight', '用法：run_demo(''project'', ''20250920'')');
    end
    out.r = start_process(varargin{1}, varargin{2:end});

% ---------------------------------------------------------------- bubble
case 'bubble'
    if isempty(varargin)
        error('run_demo:needNight', '用法：run_demo(''bubble'', ''20250920'')');
    end
    out.r = q3_bubble(varargin{1}, varargin{2:end});

% ---------------------------------------------------------------- moon
case 'moon'
    if isempty(varargin)
        error('run_demo:needNight', '用法：run_demo(''moon'', ''20260326'')');
    end
    out.r = q5_moon(varargin{1}, varargin{2:end});

% ---------------------------------------------------------------- verify
case 'verify'
    % 只跑包内自检脚本，不碰主链
    verify_paths();
    verify_paper_numbers();
    verify_checkcode();
    out.r = run_extra_selftests();

otherwise
    error('run_demo:badMode', ...
        ['未知模式 "%s"。可用：info / check / calib / phase1 / project / bubble / moon / verify'], mode);
end
end

% ==========================================================================
function show_card(tag, matFile)
%SHOW_CARD 一个 .mat 的关键参数摘要卡。
fprintf('\n  [%s]  %s\n', tag, matFile);
if exist(matFile, 'file') ~= 2
    fprintf('      （文件不存在）\n'); return
end
S = load(matFile);
if isfield(S, 'calMat'), m = S.calMat; else, m = S.cal; end
if isfield(m, 'x0')
    fprintf('      天顶像点 (x0,y0) = (%.3f, %.3f)   [0-based]\n', m.x0, m.y0);
end
if isfield(m, 'f')
    fprintf('      径向 r(z) = %.6g·(z %+.6g·z³)，z 用弧度\n', m.f, m.a);
end
if isfield(m, 'rot')
    fprintf('      方位定向角 rot   = %+.4f°\n', m.rot);
end
if isfield(m, 'rms_px')
    fprintf('      拟合残差 rms     = %.3f px   (%d 帧 / %d 夜)\n', ...
            m.rms_px, m.n_frames, m.n_nights);
end
if isfield(m, 'z_rim_deg')
    fprintf('      盘边天顶角 z_rim = %.2f°\n', m.z_rim_deg);
end
end

% ==========================================================================
function r = run_extra_selftests()
%RUN_EXTRA_SELFTESTS 跑 src/verify/ 下的其它自检（有就跑，没有就跳过）。
r = struct();
here = fileparts(mfilename('fullpath'));
d = dir(fullfile(here, 'src', 'verify', '*.m'));
for k = 1:numel(d)
    nm = d(k).name(1:end - 2);
    if any(strcmp(nm, {'verify_paths', 'verify_paper_numbers', 'verify_checkcode'})), continue; end
    fprintf('\n### %s\n', nm);
    try
        if is_script_file(fullfile(here, 'src', 'verify', [nm '.m']))
            eval(nm);              % 脚本型自检（如 moon_calib_selftest）自己会打印
            r.(nm) = 'ok';
        else
            fh = str2func(nm);     % ★ 不能写 str2func(nm)()：MATLAB 不支持对返回值就地调用
            r.(nm) = fh();
        end
    catch ME
        fprintf('  [自检失败] %s\n', ME.message);
        r.(nm) = ME.message;
    end
end
end

% ==========================================================================
function tf = is_script_file(p)
%IS_SCRIPT_FILE 该 .m 是脚本（没有顶层 function 定义）而不是函数文件。
fid = fopen(p, 'r', 'n', 'UTF-8');
if fid < 0, tf = true; return; end
cu = onCleanup(@() fclose(fid));
tf = true;
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    s = strtrim(ln);
    if isempty(s) || s(1) == '%', continue; end
    tf = ~strncmp(s, 'function', 8);
    break
end
end
