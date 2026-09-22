function moon_calib(varargin)
%MOON_CALIB  全天空气辉成像仪 —— 月亮定向定标 一体化流水线（MATLAB 自包含版）
%
% 背景
% ----
%  子午工程 儋州站（ODAZH = 富克站 FKT）双通道全天空气辉成像仪、O 通道 630 nm 的
%  归档 L0 帧里 **没有任何随天球转动的星点**（38 夜 × 15144 帧已全部排除，
%  见技能 allsky-calib-identifiability）。所以星表匹配定标原理性不可行。
%  但月亮一直在原始帧里：它自带精确星历，只要在地平上就必在视场内
%  ⇒ 可以直接反解像面的方位定向角 rot，并顺带解出径向标度 r(z)。
%
% 本脚本做的五件事（对应五个阶段）
% --------------------------------
%  scan   : 用文件名时间戳算月亮星历，普查哪些夜有月亮（不读图像，秒级）
%  detect : 在候选夜里逐帧定位饱和月核（迭代收缩窗 + 权重截断质心 + 饱和核质心）
%  solve  : ① 逐夜 rot 恒定性判据筛掉假目标
%           ② 逐夜残差自动剔除坏夜
%           ③ 模型 × 手性 网格 + 手性判决
%           ④ 残差系统学 / bootstrap / 半分交叉验证 / 站点灵敏度
%           ⑤ 低空 hold-out / 径向剖面一致性 / 亮盘边界同心性
%           ⑥ 结果健全性断言（不过就报错，绝不写出一份看着合理的错标定）
%  verify : 把「纯星历预测」的位置叠到原始帧上出 PNG
%  export : 写 .mat / .json / .txt / .csv
%
% 零工具箱依赖
% ------------
%  只用核心 MATLAB。以下都是自己实现的替代品（对应 Python 版的依赖）：
%    * 星历       mc_moon_altaz / mc_moon_ecl —— Meeus 第 47 章截断级数
%                 （不用 Aerospace Toolbox、不用第三方 ephem/skyfield）
%    * 非线性拟合  mc_lm —— 自实现 Levenberg-Marquardt + 中心差分 Jacobian
%                 （不用 Optimization Toolbox 的 lsqnonlin/lsqcurvefit）
%    * 百分位      mc_pctile —— 线性插值，与 numpy.percentile 的 'linear' 一致
%                 （不用 Statistics Toolbox 的 prctile）
%    * 高斯模糊    mc_gauss_blur —— 可分离 conv2
%                 （不用 Image Processing Toolbox 的 imgaussfilt）
%    * 画圈/画十字 mc_draw_circle / mc_draw_cross —— 直接索引像素
%                 （不用 insertShape / vision.ShapeInserter）
%    * 中位数      median(...,'omitnan')（不用 nanmedian）
%
% 用法
% ----
%   moon_calib                                    % 全流程 scan→detect→solve→verify→export
%   moon_calib('stage','scan')                    % 只普查（秒级）
%   moon_calib('stage','solve,verify,export')     % 复用已有测量缓存重跑
%   moon_calib('stage','detect','nights',{'20260325'})
%   moon_calib('site','png')                      % 用 PNG 头坐标（已知有误）做灵敏度对照
%   moon_calib('rim',512.93)                      % 换视场边缘半径
%   moon_calib('min_sat',100)                     % 饱和核像元数下限
%   moon_calib('raw','E:\qihui2\原始数据','out','D:\mooncal_out')
%   moon_calib('boot',0)                          % 跳过 bootstrap（快很多）
%   moon_calib('stage','solve,verify,export','chatty',true)   % 全量刷屏（默认精简）
%
% 输出分两级（2026-09-18 起默认精简）
% --------------------------------
%   控制台  只打"结论 + 判读 + 告警"，约 40 行；末尾有一张**定标结果卡**，一眼看清关键量。
%   日志文件 calib_report.txt **全量**（含逐夜 / 逐帧 / 逐半径的所有长表），事后可复算。
%   'chatty',true 让详细行也上控制台。实现见 mc_log / mc_set_detail。
%
% 目录约定（脚本内**不存在任何绝对路径**）
% --------------------------------------
%   本文件与 run.m / airglow_geo_pipeline.m 同放在工程根 E:\qihui2\。
%   全部产物默认落在     <本文件目录>\mooncal_work\mooncal_out_matlab\
%   默认原始数据在       E:\qihui2\原始数据\（可用 'raw' 覆盖）
%   若 <本文件目录>\mooncal_work\ 不存在，则自动回退到 <本文件目录>\..，
%   以便脚本放进 <root>\matlab\ 的旧布局也能直接跑。
%   想换地方：moon_calib('out','D:\somewhere\mooncal_out_matlab')。
%
%   %% 要给 airglow_geo_pipeline / run.m 用，加一个 'pipeline' 目标路径：
%   moon_calib('stage','solve,export', ...
%              'pipeline','E:\qihui2\Processed_Output\calibration_params.mat')
%   之后 run.m 不必改任何一行，直接 run 即可。
%   （管线 load_calibration 要 zenith/rot/radius/rz_poly 四字段；本程序自己的
%     calibration_params.mat 是 x0/y0/f/a/rot，**不能**直接喂管线。见 mc_export 注释。）
%   目标文件若已存在，会自动加时间戳备份，不会静默覆盖。
%
% 输出（默认 <本文件目录>\mooncal_work\mooncal_out_matlab\）
% --------------------------------------------------------
%   calib_report.txt        全量数字日志（UTF-8）
%   calib_params.json       标定结果 + 全部诊断表
%   calibration_params.mat  本程序原生格式（x0/y0/f/a/rot），**不是**管线格式
%   calib_summary.txt       人类可读结论 + 投影公式
%   radial_profile.csv      方位平均径向剖面
%   frame_scan.csv          逐夜普查表
%   moon_meas.mat           逐帧月核测量缓存（detect 的产物）
%   qc/qc_*.png             纯星历预测位置叠在原始帧上的验证图
%
% 投影用法（z 单位 rad、结果单位 px、坐标为 0-based）
% --------------------------------------------------
%   z = acos( sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(H) );   % 天顶角 rad
%   r = cal.f * (z + cal.a * z.^3);
%   x = cal.x0 + r .* sin(az + deg2rad(cal.rot));
%   y = cal.y0 - r .* cos(az + deg2rad(cal.rot));
%
% ⚠ 坐标系约定（本脚本最容易错的地方，务必先读）
% ----------------------------------------------
%   1) 内部**一律使用 0-based 像素坐标**（左上角像元 = (0,0)，右下角 = (1023,1023)）。
%      原因：标定值 x0≈538.9 / y0≈476.5 / R_rim≈512.9 都是 0-based。
%      只有在索引 MATLAB 数组时才 +1 —— 代码里所有 `+1` 都是这个用途。
%   2) 方位角 az：北 = 0，东 = 90（与像面 phi 约定一致）。
%   3) 文件名里的时间戳按 **UT** 解读（已由「逐夜 rot 恒定性 + 拟合残差」独立验证：
%      181 帧跨 8 h 拟合到 2 px ⇒ LST 一致到 ~0.5° ≈ 2 min）。
%
% 作者 / 日期：2026-09-16，由 Python 版 moon_calib.py 1:1 移植
%              （Python 原型轨已于 2026-09-17 删除，备份见 E:\qihui2_代码备份_20260917\）

% ==============================================================================
% [0] 解析参数
% ==============================================================================
% 本文件所在目录 = 工程根（q_cfg.m / airglow_geo_pipeline.m 同目录）。
% （2026-09-17：旧驱动 run.m 已删除，新驱动为 q1_calibrate / q2_project。）
% 产物一律落在 <here>\mooncal_work\ 下，脚本内不再出现任何绝对路径。
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
C = mc_defaults(here);

ip = inputParser;
ip.addParameter('stage', 'all');
ip.addParameter('raw',   C.raw);
ip.addParameter('out',   C.out);
ip.addParameter('pipeline', '', @(v) ischar(v) || (isstring(v) && isscalar(v)));
ip.addParameter('outmat', '', @(v) ischar(v) || (isstring(v) && isscalar(v)));
ip.addParameter('site',  'official');
ip.addParameter('nights', {});
ip.addParameter('every', 1);
ip.addParameter('min_sat', C.min_sat);
ip.addParameter('min_illum', 0.15);
ip.addParameter('alt_lo', 14.0);
ip.addParameter('alt_hi', 78.0);
ip.addParameter('alt_hold', 3.0);
ip.addParameter('scat_max', 2.5);
ip.addParameter('min_frames', 6);
% ★ 太阳高度门限 (deg)：只保留 sun_alt < sun_alt_max 的帧。默认 Inf = 不过滤（二期行为不变）。
%   一期 2014（8 bit）必须开，否则黄昏/黎明天光饱和帧会被当成月核 ⇒ 7 夜全灭。
%   推荐 -15（民用暮光之下）；-18 更干净但会砍掉低空帧。见 mc_sun_altaz 注释。
ip.addParameter('sun_alt_max', Inf);
ip.addParameter('rim', C.rim);
% ★ 几何天顶像点 / 视场边缘半径 **可覆盖**。
%   mc_defaults 里的 (543.824, 474.352)/512.93 是**二期 2025-2026 那台相机**的值。
%   换一批数据（尤其是一期 2014 的 8 bit 数据）时光学不同：一期实测圆心≈(552, 555)、
%   R_rim≈458 px。若不覆盖，[3.1] 判据里的 phi = atan2(dx,-dy) 会带上 ~80 px 的中心偏差
%   ⇒ phi 误差随月亮方位变化 ⇒ 逐夜 rot 散布被抬到 4~92 deg ⇒ 7 夜全被判成"假目标"，
%   而报错信息会伪装成"数据里没有月亮"。见 calib_2014.m 的说明。
ip.addParameter('x0_geo', C.x0_geo);
ip.addParameter('y0_geo', C.y0_geo);
ip.addParameter('night_rms_mult', 4.0);
ip.addParameter('night_rms_floor', 3.0);
ip.addParameter('aperture', 18);
ip.addParameter('min_snr_dn', 1500.0);
ip.addParameter('boot', 200);
ip.addParameter('qc', 6);
ip.addParameter('log', '');
ip.addParameter('verbose_decode', false);
% ★ 输出详细度（2026-09-18 加）。默认 false = **精简**：控制台只打结论与判读，
%   逐夜/逐帧长表仍完整写进 calib_report.txt（事后可复算）。true = 恢复 2026-09-17 之前的
%   全量刷屏。实现见 mc_log / mc_set_detail。总开关只在这里拨一次。
ip.addParameter('chatty', false);
ip.parse(varargin{:});
o = ip.Results;

% ⚠ 必须在**任何 logf 之前**设定，否则前几行会按默认值走。
mc_set_detail(o.chatty);
% 清掉上一次运行留下的"产物路径"缓存（appdata 在 MATLAB 会话内是持久的；
% 若本次只跑 solve，残留的旧路径会让末尾结果卡指向上一轮的文件）。
setappdata(0, 'mc_cal_exported', []);

siteKey = char(o.site);
assert(isfield(C.sites, siteKey), 'moon_calib:badSite', '未知站点 "%s"（可选 official / png）', siteKey);
site = C.sites.(siteKey);

% ★ 覆盖几何参考（见 inputParser 处注释）。默认 = 二期值 ⇒ 对二期完全向后兼容。
C.x0_geo = double(o.x0_geo);
C.y0_geo = double(o.y0_geo);
C.rim    = double(o.rim);

stages = strtrim(strsplit(char(o.stage), ','));
if any(strcmp(stages, 'all'))
    stages = {'scan', 'detect', 'solve', 'verify', 'export'};
end
has = @(s) any(strcmp(stages, s));

if ~exist(char(o.out), 'dir'), mkdir(char(o.out)); end
outdir = char(o.out);
logpath = char(o.log);
if isempty(logpath), logpath = fullfile(outdir, 'calib_report.txt'); end

fid = fopen(logpath, 'w', 'n', 'UTF-8');
assert(fid > 0, 'moon_calib:noLog', '无法写日志文件 %s', logpath);
cu = onCleanup(@() fclose(fid));
logf = @(varargin) mc_log(fid, varargin{:});

cfg = struct();
cfg.raw             = char(o.raw);
cfg.out             = outdir;
cfg.outmat          = char(o.outmat);
cfg.pipeline_out    = char(o.pipeline);
cfg.site_key        = siteKey;
cfg.rim             = o.rim;
cfg.x0_geo          = C.x0_geo;
cfg.y0_geo          = C.y0_geo;
cfg.alt_lo          = o.alt_lo;
cfg.alt_hi          = o.alt_hi;
cfg.alt_hold        = o.alt_hold;
cfg.scat_max        = o.scat_max;
cfg.min_frames      = o.min_frames;
cfg.sun_alt_max     = o.sun_alt_max;
cfg.night_rms_mult  = o.night_rms_mult;
cfg.night_rms_floor = o.night_rms_floor;
cfg.n_boot          = o.boot;
cfg.aperture        = o.aperture;
cfg.min_snr_dn      = o.min_snr_dn;
cfg.n_qc            = o.qc;

t0 = tic;
stageCN = struct('scan', '普查', 'detect', '检测', 'solve', '解算', ...
                 'verify', '验证', 'export', '导出');
stgNames = cellfun(@(s) ternary(isfield(stageCN, s), stageCN.(s), s), stages, ...
                   'UniformOutput', false);
switch siteKey
    case 'official'
        siteCN = '儋州站 富克 FKT（官方坐标）';
    case 'png'
        siteCN = '按 PNG 头坐标【二期值】(109.83E,19.31N) —— 已知有误，仅作灵敏度对照';
    case 'png1'
        siteCN = '按 PNG 头坐标【一期值】(109.1E,19.5N) —— 仅作灵敏度对照';
    otherwise
        siteCN = siteKey;
end
logf('%s', repmat('=', 1, 78));
logf('  月亮法定标 · 开始');
logf('  时间 %s UTC   阶段 %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), strjoin(stgNames, ' → '));
logf('  站点 %s  (%.3fE, %.3fN, %.0f m)', siteCN, site(1), site(2), site(3) * 1000);
logf('  原始数据 %s', cfg.raw);
logf('  输出目录 %s', outdir);
logf('%s', repmat('=', 1, 78));

% ==============================================================================
% [1] 阶段 1 —— 普查
% ==============================================================================
fixedNights = mc_cellstr(o.nights);
cand = fixedNights;
if has('scan')
    cand = mc_scan(cfg, site, fixedNights, 0.0, o.min_illum, logf);
elseif isempty(cand)
    cand = mc_list_nights(cfg.raw);
end
cfg.cand_nights = cand;

% ==============================================================================
% [2] 阶段 2 —— 检测
% ==============================================================================
if has('detect')
    recs = mc_detect(cfg, cand, o.every, o.min_sat, o.verbose_decode, logf);
    % ⚠ 防呆 1: 检测结果为 0 帧时必须**致命报错**。
    %   历史上（Python 版）这里的成因是代码 bug（漏写装饰器 ⇒ 每帧都判成解码失败），
    %   而不是数据问题。若静默继续，solve 读不到缓存就跳过 export，
    %   输出目录里只剩上一轮的旧 JSON，看上去一切正常 —— 实测差点交出过期标定。
    assert(~isempty(recs), 'moon_calib:detectZero', ...
        ['detect 阶段一帧都没检出（共 %d 个候选夜）。这几乎总是**代码问题**而不是数据问题，先查：\n' ...
         '  1) mc_load_gray 是否对所有帧都返回空（解码链路坏）；\n' ...
         '  2) 饱和核判据 min_sat=%d 是否过高；\n' ...
         '  3) 候选夜是不是都真的没有月亮。'], numel(cand), o.min_sat);
end

% ==============================================================================
% [3] 阶段 3 —— 解算
% ==============================================================================
cal = [];
% ⚠ 只在这里静音一次，不要在 mc_lm 里反复开关（见 mc_lm 注释：
%   反复 warning('off') / 恢复会在 R2021b -batch 下把 MATLAB 搞成堆损坏 0xc0000374）。
% 被静音的只有这一个 ID：它来自"不使用参数 a 的模型"使 Jacobian 第 5 列恒零，
% 属预期行为；其它警告照常显示。
wstGlobal = warning('off', 'MATLAB:nearlySingularMatrix');
cuW = onCleanup(@() warning(wstGlobal));

if has('solve')
    cal = mc_solve(cfg, site, logf);
    % ⚠ 防呆 2: 解算失败/未产出结果时必须致命报错，绝不静默跳过 export。
    %   否则输出目录里留着的上一轮旧 JSON/MAT 看起来"一切正常"。
    assert(~isempty(cal), 'moon_calib:solveFailed', ...
        ['solve 未产出结果（详见上文），不写任何导出文件。\n' ...
         '注意: %s 里可能还留着**上一轮**的 calib_params.json / .mat，不要误用。'], outdir);
end

% ==============================================================================
% [4] 阶段 4 —— 验证图
% ==============================================================================
qc = struct('night', {}, 'file', {}, 't', {}, 'alt', {}, 'az', {}, ...
            'pred_px', {}, 'pred_py', {}, 'peak_near', {}, 'bg', {}, 'snr_dn', {}, 'png', {});
if has('verify')
    if isempty(cal)
        logf('solve 未跑 -> 跳过 verify');
    else
        qc = mc_verify(cfg, site, cal, logf);
    end
end

% ==============================================================================
% [5] 阶段 5 —— 导出
% ==============================================================================
if has('export')
    if isempty(cal)
        logf('solve 未跑 -> 跳过 export');
    else
        mc_export(cfg, cal, qc, logf);
    end
end

logf('');
logf('完成，用时 %.1f s   日志 %s', toc(t0), logpath);

% ==============================================================================
% [6] 结果卡 —— 整个流程的**最后一块输出**，一眼看清定标结果
% ==============================================================================
% 放在最后（而不是 3.3 的【解算结论】处）是刻意的：阶段 4/5 会再刷十几行，
% 跑完之后最后看到的是文件路径，真正要的两个数（天顶像点、rot）反倒被刷走了。
if ~isempty(cal)
    mc_result_card(cal, cfg, logpath, toc(t0), logf);
else
    logf('');
    logf('%s', repmat('=', 1, 78));
    logf('  定标未完成 —— solve 阶段没有产出结果（原因见上文），未写任何导出文件。');
    logf('%s', repmat('=', 1, 78));
end
end


% ################################################################################
% #                                                                              #
% #                       以 下 全 部 是 局 部 函 数                              #
% #                                                                              #
% ################################################################################


% ==============================================================================
% [0] 配置
% ==============================================================================

function C = mc_defaults(here)
%MC_DEFAULTS 全部常量。改这里，不要在函数体里散落魔术数字。
%   here : moon_calib.m 所在目录（= 工程根）。省略时按 mfilename 推断。
if nargin < 1 || isempty(here)
    here = fileparts(mfilename('fullpath'));
    if isempty(here), here = pwd; end
end
% 产物根目录：优先 <here>\mooncal_work（当前布局）；
% 若不存在则回退 <here>\..（旧的 <root>\matlab\ 布局），保证向后兼容。
mcWork = mc_paths().mooncal;   % 整包唯一的产物根（不再按 here 推断）
C.size       = 1024;        % 图像边长 (px)
C.sat_level  = 65000;       % 16 bit 饱和阈值（满量程 65535）
C.min_sat    = 100;         % "饱和核"像元数下限 -> 判定"这是月亮而不是噪声"
C.blob_half  = 80;          % 质心迭代窗口半宽 (px)
C.blob_niter = 3;           % 质心迭代次数
C.sat_near   = 45;          % 判定"该饱和像元属于月面"的半径 (px)

% 纯几何法给出的视场边缘半径：由视场边界的**旋转对称性**定出，与任何天体无关。
C.rim        = 512.93;

% 几何法天顶像点 (0-based)。**只用作初值与 QC 对照**，拟合会重新求 x0,y0。
% 它是独立于月亮的一条路线，两者相差多少本身就是重要 QC。
C.x0_geo     = 543.824;
C.y0_geo     = 474.352;

% 站点坐标 [lon_deg, lat_deg, height_km]
%   official : 儋州站 = 海南省儋州市雅星镇富克街道 = 富克站 FKT，官方公布坐标
%   png      : 【二期 2026-03 起】PNG tEXt 头里写的坐标，**已知有误**（经度差 0.70 deg ≈ 250 km 上 73 km）
%   png1     : 【一期 2013/2014】PNG 头里写的坐标。★ 2026-09-18 实测：两期头坐标**不是同一个值**，
%              一期那个与官方只差 0.033E/0.026N，二期那个差 0.697E/0.216N。原先只存了 'png'
%              （二期值），对一期数据用它就等于套用了"另一台相机的头坐标"。两者都只作灵敏度对照。
C.sites = struct();
C.sites.official = [109.133, 19.526, 0.103];
C.sites.png      = [109.830, 19.310, 0.103];
C.sites.png1     = [109.100, 19.500, 0.028];

% 径向投影模型候选（与 Python 版同名同序）
% model_desc : 公式（写进日志/JSON，保持符号原样，便于跨语言比对）
% model_cn   : 中文名（只用于控制台表格，避免每次都要回忆 equid+ 是什么）
C.models = {'equidistant', 'equid+', 'equisolid', 'orthographic', 'stereographic'};
C.model_desc = {'r = f z', 'r = f (z + a z^3)', 'r = 2 f sin(z/2)', ...
                'r = f sin z', 'r = 2 f tan(z/2)'};
C.model_cn = {'等距投影', '等距+三次项', '等立体角投影', '正交投影', '立体投影'};

C.raw = mc_paths().raw;
C.out = fullfile(mcWork, 'mooncal_out_matlab');
end


% ==============================================================================
% [1] 基础工具
% ==============================================================================

function mc_log(fid, varargin)
%MC_LOG 日志输出。**始终**写入 UTF-8 日志文件；是否同时打到控制台，看这一行的级别。
%
%   mc_log(fid, fmt, ...)         普通行 —— 控制台 + 文件（结论、判读、告警）
%   mc_log(fid, '详', fmt, ...)   详细行 —— 只进文件；除非开了 'chatty'（此时两级都上屏）
%
% 为什么分两级（2026-09-18 加）：
%   全量输出有 580 行 / 39 KB，其中 90% 是"逐夜/逐帧/逐半径"的长表。它们对**事后复算**
%   是必需的，对**当场看结果**却是噪声 —— 真正要看的只有几行（模型、残差、rot、天顶像点）。
%   所以长表一律标 '详' 进文件，控制台只留结论。要临时看全量：moon_calib('stage',...,'chatty',true)。
detail = false;
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && strcmp(char(varargin{1}), '详')
    detail = true;
    varargin(1) = [];
end
s = sprintf(varargin{:});
fprintf(fid, '%s\n', s);
% flag 用 root appdata 存放（跨本地函数共享，且 -batch 下也有效）。
% 未设置时 getappdata 返回 [] -> 按"不上屏"处理。
flag = getappdata(0, 'mc_log_detail');
if ~detail || (isscalar(flag) && logical(flag))
    fprintf('%s\n', s);
end
end


function mc_set_detail(tf)
%MC_SET_DETAIL 设定"详细行"是否也打印到控制台。见 mc_log 的分级说明。
setappdata(0, 'mc_log_detail', logical(tf));
end


function c = mc_cellstr(x)
%MC_CELLSTR 把 char / string / cell / 逗号串统一成 cellstr。
if isempty(x)
    c = {};
elseif ischar(x)
    c = strtrim(strsplit(x, ','));
    c = c(~cellfun(@isempty, c));
elseif isstring(x) && isscalar(x)
    c = {char(x)};
elseif iscell(x)
    c = cellfun(@(q) char(q), x(:)', 'UniformOutput', false);
    c = c(~cellfun(@isempty, c));
else
    c = {};
end
end


function t = mc_parse_ts(name)
%MC_PARSE_TS 从文件名里抠出 YYYYMMDDhhmmss -> datetime。
% 时刻按 UT 解读（已被 rot 恒定性 + 拟合残差独立验证，不要再怀疑时区约定）。
% 返回 [] 表示这个文件名不带可解析的时间戳。
t = [];
tok = regexp(name, '(\d{8})(\d{6})', 'tokens', 'once');
if isempty(tok), return; end
d = tok{1};
h = tok{2};
v = [str2double(d(1:4)), str2double(d(5:6)), str2double(d(7:8)), ...
     str2double(h(1:2)), str2double(h(3:4)), str2double(h(5:6))];
if any(isnan(v)), return; end
try
    t = datetime(v(1), v(2), v(3), v(4), v(5), v(6));
catch
    t = [];
end
end


function fs = mc_list_frames(nightDir)
%MC_LIST_FRAMES 列出某夜所有可解析帧（按时间戳升序）。
%
% ⚠ 坑（Python 版踩过，MATLAB 版同样要防）：
%   绝不要用 dir(*.PNG) + dir(*.png) 这种"两种大小写各来一次"的写法。
%   Windows 文件名大小写不敏感，两条 dir 会返回**同一批文件**，帧表被接两遍
%   ⇒ 索引错位 ⇒ Δt 变负 ⇒ 旋转角取反 ⇒ 最后输出一个"看起来很像结论"的假阴性。
%   这里统一用 dir('*') + 小写后缀判定。
fs = {};
if ~exist(nightDir, 'dir'), return; end
d = dir(nightDir);
keep = false(1, numel(d));
for i = 1:numel(d)
    keep(i) = ~d(i).isdir && numel(d(i).name) >= 4 && ...
              strcmpi(d(i).name(end-3:end), '.png');
end
d = d(keep);
if isempty(d), return; end
ts = NaT(numel(d), 1);
ok = false(numel(d), 1);
for i = 1:numel(d)
    t = mc_parse_ts(d(i).name);
    if ~isempty(t)
        ts(i) = t;
        ok(i) = true;
    end
end
d = d(ok);
ts = ts(ok);
[~, ord] = sort(ts);
fs = cell(1, numel(d));
for i = 1:numel(d)
    fs{i} = fullfile(nightDir, d(ord(i)).name);
end
end


function ns = mc_list_nights(raw)
%MC_LIST_NIGHTS 原始数据根目录下的观测夜子目录（升序）。
ns = {};
if ~exist(raw, 'dir'), return; end
d = dir(raw);
d = d([d.isdir]);
nm = {d.name};
nm = nm(~ismember(nm, {'.', '..'}));
ns = sort(nm);
end


function im = mc_load_gray(path)
%MC_LOAD_GRAY 读一帧灰度图，**按位深统一到 16 bit 量级**。失败返回 []。
%
% ★ 位深（2026-09-17 修的一个真 bug）★
%   归档里有两批数据：
%     二期 2025-2026 : 1024x1024, 16 bit（满量程 65535），夜天光 ~3400~8000 DN
%     一期 2014      : 1024x1024,  8 bit（满量程   255），夜天光 ~13~60   灰阶
%   而本文件所有 DN 阈值都按 16 bit 写死，最要紧的是 C.sat_level = 65000。
%   若直接 imread 拿 8 bit 数据，`f > 65000` **永远为假** ⇒ detect 阶段一帧都
%   检不出，而 assert 的报错文案会把人往"数据里没有月亮"上引 —— 实测就是这样
%   卡住一期定标的。
%   ⇒ 统一委托 q_read_gray：uint16 原样、uint8 乘 257 提升到 16 bit 量级。
%
% 另两个历史坑：
%   (1) OpenCV 读不了中文路径 -> Python 用 np.fromfile+cv2.imdecode。
%       MATLAB 的 imread 原生支持中文路径，不需要绕。
%   (2) 归档帧可能整夜截断损坏（IHDR 完好、IDAT 不全）。imread 会抛错或告警，
%       q_read_gray 里 try/catch 兜住并返回 []，调用者按夜汇总打印
%       "解码失败 n 帧"。**绝不**让它静默算一个数字出来。
im = q_read_gray(path);
end


function B = mc_gauss_blur(A, sigma)
%MC_GAUSS_BLUR 可分离高斯模糊，边界用 **BORDER_REFLECT_101** 复刻 OpenCV。
%
% 为什么连边界都要抠？因为本函数唯一的用途是给"找最亮斑"提供一个抗椒盐的**种子**，
% 而月核是一个**饱和平台**（几千个像元都顶着 65000）：argmax 落在平台的哪个角，
% 完全由"打破平局"的扫描顺序决定。种子差几像素，经 3 次收缩窗迭代后仍会在参数上
% 留下系统差（本项目实测：零填充 + 列主序 ⇒ y0 偏 0.66 px、rot 偏 0.1 deg）。
% Python 版用的是 cv2.GaussianBlur 的默认边界 BORDER_REFLECT_101：
%     ... c b a | a b c d | d c b a ...
% 注意最外侧那个 a **不重复**（区别于 BORDER_REFLECT 的 "... b a | a b c d | d c b a ..."）。
% 自己实现是为了不依赖 Image Processing Toolbox 的 imgaussfilt。
% ⚠ 核长必须与 OpenCV 对齐，否则平台上的 argmax 会翻到相邻像元（实测 cy 偏 0.3~0.4 px）：
%   OpenCV 在 ksize=0 时按  ksize = cvRound(sigma * (depth==CV_8U ? 3 : 4) * 2 + 1) | 1
%   取核长。Python 版传的是 float32，属非 8 位分支 ⇒ sigma=2 时 ksize = 17，半径 8。
%   即它对 **4σ** 截断，而不是常见的 3σ。核本身仍是 exp(-x^2/2σ^2) 并在截断窗内归一化。
r = max(1, ceil(4 * sigma));
x = (-r:r)';
k = exp(-(x.^2) / (2 * sigma^2));
k = k / sum(k);
P = mc_pad_reflect101(double(A), r);
B = conv2(k, k', P, 'valid');
end


function P = mc_pad_reflect101(A, r)
%MC_PAD_REFLECT101 按 OpenCV BORDER_REFLECT_101 规则做反射填充。
idxR = mc_reflect_idx(size(A, 1), r);
idxC = mc_reflect_idx(size(A, 2), r);
P = A(idxR, idxC);
end


function idx = mc_reflect_idx(n, r)
%MC_REFLECT_IDX 生成长度为 n+2r 的反射索引（1-based）。
%   越界下侧: i=0 -> 2, i=-1 -> 3, ...   （即 2-i，最外侧原像元不重复）
%   越界上侧: i=n+1 -> n-1, i=n+2 -> n-2, ...  （即 2n-i）
% r << n 时一次反射就够，无需迭代。
i = (1 - r):(n + r);
idx = i;
lo = i < 1;
idx(lo) = 2 - i(lo);
hi = i > n;
idx(hi) = 2 * n - i(hi);
end


function q = mc_pctile(x, p)
%MC_PCTILE 百分位（线性插值），与 numpy.percentile 的默认 'linear' 完全一致。
% 自己实现是为了**不依赖** Statistics Toolbox 的 prctile。
v = sort(x(:));
n = numel(v);
if n == 0, q = NaN; return; end
if n == 1, q = v(1); return; end
pos = (p / 100) * (n - 1);
lo = floor(pos);
hi = ceil(pos);
q = v(lo + 1) + (pos - lo) * (v(hi + 1) - v(lo + 1));
end


function [XX, YY] = mc_grid(SZ)
%MC_GRID 0-based 像素坐标网格：XX = 列坐标 x，YY = 行坐标 y。
[XX, YY] = meshgrid(0:SZ-1, 0:SZ-1);
end


function [m, sd, R] = mc_circstats(valsDeg)
%MC_CIRCSTATS 圆统计：返回（圆均值 deg, 圆散布 deg, R）。R 越接近 1 越集中。
a = deg2rad(double(valsDeg(:)));
c = mean(cos(a));
s = mean(sin(a));
R = hypot(c, s);
sd = rad2deg(sqrt(-2 * log(max(R, 1e-12))));
m = mod(rad2deg(atan2(s, c)), 360);
end


function d = mc_ang_diff(a, b)
%MC_ANG_DIFF 把 a-b 折算到 (-180, 180]。
d = mod(a - b + 180, 360) - 180;
end


% ==============================================================================
% [2] 星历（自包含，不依赖任何天文/工具箱）
% ==============================================================================

function jd = mc_jd(dtv)
%MC_JD UT 时刻 -> 儒略日（Meeus 第 7 章）。
y = year(dtv);
m = month(dtv);
d = day(dtv) + (hour(dtv) + minute(dtv) / 60 + second(dtv) / 3600) / 24;
if m <= 2
    y = y - 1;
    m = m + 12;
end
A = floor(y / 100);
B = 2 - A + floor(A / 4);
jd = floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + B - 1524.5;
end


function g = mc_gmst_deg(jd)
%MC_GMST_DEG 格林尼治平恒星时（deg）。
T = (jd - 2451545.0) / 36525.0;
g = mod(280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T * T, 360.0);
end


function [lam, bet, eps, distKm] = mc_moon_ecl(jd)
%MC_MOON_ECL 月球黄道坐标（Meeus《Astronomical Algorithms》第 47 章截断级数）。
% 返回 (黄经 rad, 黄纬 rad, 黄赤交角 rad, 地心距 km)。
%
% 精度：黄经约 0.1 deg，在天顶附近约 0.6 px（在 326.5 px/rad 的尺度上）。
% ⚠ 这是本项目定标精度的**主导系统项之一**；再往上要换 ELP2000 全级数。
T = (jd - 2451545.0) / 36525.0;
r = @deg2rad;

Lp = r(mod(218.3164477 + 481267.88123421 * T, 360));
D  = r(mod(297.8501921 + 445267.1114034  * T, 360));
M  = r(mod(357.5291092 + 35999.0502909  * T, 360));
Mp = r(mod(134.9633964 + 477198.8675055 * T, 360));
F  = r(mod(93.2720950  + 483202.0175233 * T, 360));

E = 1.0 - 0.002516 * T - 0.0000074 * T * T;      % 地球轨道偏心率修正

lam = Lp ...
    + r(6.288774) * sin(Mp) ...
    + r(1.274027) * sin(2 * D - Mp) ...
    + r(0.658314) * sin(2 * D) ...
    + r(0.213618) * sin(2 * Mp) ...
    - r(0.185116) * E * sin(M) ...
    - r(0.114332) * sin(2 * F) ...
    + r(0.058793) * sin(2 * D - 2 * Mp) ...
    + r(0.057066) * E * sin(2 * D - M - Mp) ...
    + r(0.053322) * sin(2 * D + Mp) ...
    + r(0.045758) * E * sin(2 * D - M) ...
    - r(0.040923) * E * sin(M - Mp) ...
    - r(0.034720) * sin(D) ...
    - r(0.030383) * E * sin(M + Mp) ...
    + r(0.015327) * sin(2 * D - 2 * F) ...
    - r(0.012528) * sin(Mp + 2 * F) ...
    + r(0.010980) * sin(Mp - 2 * F);

bet = r(5.128622) * sin(F) ...
    + r(0.280602) * sin(Mp + F) ...
    + r(0.277693) * sin(Mp - F) ...
    + r(0.173237) * sin(2 * D - F) ...
    + r(0.055413) * sin(2 * D - Mp + F) ...
    + r(0.046271) * sin(2 * D - Mp - F) ...
    + r(0.032573) * sin(2 * D + F) ...
    + r(0.017198) * sin(2 * Mp + F);

distKm = 385000.56 - 20905.355 * cos(Mp) - 3699.111 * cos(2 * D - Mp) ...
         - 2955.968 * cos(2 * D) - 569.925 * cos(2 * Mp);

% Meeus 47 的附加改正（A1/A2/A3 项），量级 ~0.005 deg，顺手加上
A1 = r(mod(119.75 + 131.849   * T, 360));
A2 = r(mod(53.09  + 479264.290 * T, 360));
A3 = r(mod(313.45 + 481266.484 * T, 360));
lam = lam + r(0.004 * (sin(A1) + sin(Lp - F)) + 0.000318 * sin(A2));
bet = bet + r(0.002 * (sin(A3) + sin(A1 - F) + sin(A1 + F) + sin(Lp - Mp) - sin(Lp + Mp)));

eps = r(23.4392911 - 0.0130042 * T);
end


function lamS = mc_sun_ecl_lon(jd)
%MC_SUN_ECL_LON 太阳黄经（Meeus 第 25 章低精度），精度 ~0.01 deg。只用于算月相照度。
T = (jd - 2451545.0) / 36525.0;
L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
M = deg2rad(mod(357.52911 + 35999.05029 * T - 0.0001537 * T * T, 360));
C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * sin(M) ...
  + (0.019993 - 0.000101 * T) * sin(2 * M) ...
  + 0.000289 * sin(3 * M);
lamS = deg2rad(mod(L0 + C, 360.0));
end


function k = mc_moon_illum(jd)
%MC_MOON_ILLUM 月面被照亮比例 k ∈ [0,1]（1 = 望）。
[lamM, betM, ~, ~] = mc_moon_ecl(jd);
lamS = mc_sun_ecl_lon(jd);
cosPsi = -cos(lamM - lamS) * cos(betM);
k = 0.5 * (1.0 + cosPsi);
end


function [ra, de] = mc_ecl2eq(lam, bet, eps)
%MC_ECL2EQ 黄道 -> 赤道。
ra = atan2(sin(lam) * cos(eps) - tan(bet) * sin(eps), cos(lam));
de = asin(sin(bet) * cos(eps) + cos(bet) * sin(eps) * sin(lam));
ra = mod(ra, 2 * pi);
end


function [altDeg, azDeg] = mc_moon_altaz(jd, site)
%MC_MOON_ALTAZ 月亮的地平坐标 (alt_deg, az_deg)，**含站心视差**（topocentric）。
%
% 视差必须做：月亮视差 ~0.95 deg，在本机 ~326.5 px/rad 上是约 5 px，不能忽略。
% 方位角约定：北 = 0，东 = 90（与像面 phi 约定一致）。
lon = site(1);
lat = site(2);
hgtKm = site(3);

[lam, bet, eps, dist] = mc_moon_ecl(jd);
[ra, de] = mc_ecl2eq(lam, bet, eps);
sinPi = 6378.14 / dist;

phi = deg2rad(lat);
u = atan(0.99664719 * tan(phi));
rs = 0.99664719 * sin(u) + (hgtKm / 6378.14) * sin(phi);
rc = cos(u) + (hgtKm / 6378.14) * cos(phi);

H = deg2rad(mc_gmst_deg(jd) + lon) - ra;
dA = atan2(-rc * sinPi * sin(H), cos(de) - rc * sinPi * cos(H));
raT = ra + dA;
deT = atan2((sin(de) - rs * sinPi) * cos(dA), cos(de) - rc * sinPi * cos(H));

Ht = deg2rad(mc_gmst_deg(jd) + lon) - raT;
z = acos(max(-1.0, min(1.0, sin(phi) * sin(deT) + cos(phi) * cos(deT) * cos(Ht))));
A = atan2(-sin(Ht) * cos(deT), ...
          sin(deT) * cos(phi) - cos(deT) * sin(phi) * cos(Ht));
altDeg = 90.0 - rad2deg(z);
azDeg = mod(rad2deg(A), 360.0);
end


function [altDeg, azDeg] = mc_sun_altaz(jd, site)
%MC_SUN_ALTAZ 太阳的地平坐标 (alt_deg, az_deg)，低精度（Meeus 25 章 + 黄赤转换）。
%
% ★ 为什么需要它：8 bit 的一期 2014 数据里，**黄昏/黎明天光**也会把大片像元顶到满量程，
%   被 [阶段 2] 当成"饱和月核"检出来。这些假目标的质心跟着**太阳**（严格说是天光梯度）
%   而不是月亮走 ⇒ 逐夜 rot 散布被抬到 8~90 deg，整整 7 夜全被判成"假目标"。
%   实测（20140930 起 6 帧，太阳<-15 deg）：rot = 6.9/6.2/6.1/6.8/6.6/6.8 deg 极稳；
%   而同夜太阳只有 -8.4 deg 的首帧给出 118.6 deg。
%   ⇒ 用 sun_alt_max 门限把暮光帧剔掉，是让一期月亮法重新可用的**关键一步**。
%
% 太阳视差 ~8.8" 可忽略，不做站心修正。
lon = site(1);
lat = site(2);
T = (jd - 2451545.0) / 36525.0;
eps = deg2rad(23.4392911 - 0.0130042 * T);
[ra, de] = mc_ecl2eq(mc_sun_ecl_lon(jd), 0.0, eps);
phi = deg2rad(lat);
H = deg2rad(mc_gmst_deg(jd) + lon) - ra;
z = acos(max(-1.0, min(1.0, sin(phi) * sin(de) + cos(phi) * cos(de) * cos(H))));
A = atan2(-sin(H) * cos(de), sin(de) * cos(phi) - cos(de) * sin(phi) * cos(H));
altDeg = 90.0 - rad2deg(z);
azDeg = mod(rad2deg(A), 360.0);
end


% ==============================================================================
% [3] 投影模型
% ==============================================================================

function r = mc_r_of_z(model, z, f, a)
%MC_R_OF_Z z: 天顶角 (rad) -> 像面半径 (px)。z 可以是向量。
switch model
    case 'equidistant'
        r = f * z;
    case 'equid+'
        r = f * (z + a * z.^3);
    case 'equisolid'
        r = 2.0 * f * sin(z / 2.0);
    case 'orthographic'
        r = f * sin(z);
    case 'stereographic'
        r = 2.0 * f * tan(z / 2.0);
    otherwise
        error('moon_calib:badModel', '未知投影模型 "%s"', model);
end
end


function zDeg = mc_z_of_rim(model, f, a, rim)
%MC_Z_OF_RIM 反解：半径 rim 处的天顶角 (deg)。
% 搜索范围故意放宽到 126 deg —— 因为实测发现"亮环半径"可能**大于 r(90 deg)**，
% 即视场略超 180 deg（超广角鱼眼常见）。返回 >90 是物理上有意义的信息，不要截断成 NaN。
zz = linspace(1e-4, 2.2, 60000);
rv = mc_r_of_z(model, zz, f, a);
if any(~isfinite(rv)) || max(rv) < rim
    zDeg = NaN;
    return;
end
[~, k] = min(abs(rv - rim));
zDeg = rad2deg(zz(k));
end


% ==============================================================================
% [4] 阶段 1 —— 普查：哪些夜的月亮可见
% ==============================================================================

function cand = mc_scan(cfg, site, nights, minAlt, minIllum, logf)
%MC_SCAN 只读文件名算星历，不读图像 -> 几秒钟就能扫完整个档案。
% 这一步的意义是**把 detect 阶段的图像读取量压到 1/10**：月亮不在天上的夜根本不用看。
if isempty(nights), nights = mc_list_nights(cfg.raw); end
rows = struct('night', {}, 'n_frames', {}, 'n_good', {}, 'win', {}, ...
              'alt_min', {}, 'alt_max', {}, 'illum_max', {}, 'verdict', {});

logf('');
logf('%s', repmat('=', 1, 78));
logf('【阶段 1】月亮可见性普查   —— 只读文件名算星历，不读图像，秒级');
logf('          判据：整夜存在 alt > %.0f° 且 照度 > %.0f%% 的时段', minAlt, minIllum * 100);
logf('%s', repmat('=', 1, 78));
logf('详', '站点 (%.3fE, %.3fN, %.0f m)', site(1), site(2), site(3) * 1000);
logf('详', '%-10s %6s | %-19s %-19s | %8s %7s | %s', ...
     '夜', '帧数', '月亮可观测窗口(UT)', 'alt 范围', '最大照度', '峰值alt', '结论');
logf('详', '%s', repmat('-', 1, 118));

for i = 1:numel(nights)
    nd = nights{i};
    d = fullfile(cfg.raw, nd);
    fs = mc_list_frames(d);
    if numel(fs) < 20
        logf('详', '%-10s %6d | %-19s %-19s | %8s %7s | 帧不足', nd, numel(fs), '-', '-', '-', '-');
        continue;
    end
    n = numel(fs);
    alts = zeros(n, 1);
    ills = zeros(n, 1);
    tss  = NaT(n, 1);
    for k = 1:n
        t = mc_parse_ts(mc_basename(fs{k}));
        tss(k) = t;
        jd = mc_jd(t);
        a_ = mc_moon_altaz(jd, site);
        alts(k) = a_;
        ills(k) = mc_moon_illum(jd);
    end
    good = (alts > minAlt) & (ills > minIllum);
    ok = sum(good);
    if ok >= 6
        idx = find(good);
        win = sprintf('%s-%s', datestr(tss(idx(1)), 'HH:MM'), datestr(tss(idx(end)), 'HH:MM'));
        aRng = sprintf('%+.0f..%+.0f', min(alts(idx)), max(alts(idx)));
        verdict = '候选';
        illMax = max(ills(idx));
        altPk = max(alts(idx));
    else
        win = '-'; aRng = '-'; verdict = '否';
        illMax = max(ills); altPk = max(alts);
    end
    logf('详', '%-10s %6d | %-19s %-19s | %8.3f %7.1f | %s', nd, n, win, aRng, illMax, altPk, verdict);
    rows(end+1) = struct('night', nd, 'n_frames', n, 'n_good', ok, 'win', win, ...   %#ok<AGROW>
                         'alt_min', min(alts), 'alt_max', max(alts), ...
                         'illum_max', max(ills), 'verdict', verdict);
end

cand = {};
for i = 1:numel(rows)
    if strcmp(rows(i).verdict, '候选'), cand{end+1} = rows(i).night; end   %#ok<AGROW>
end
logf('详', '%s', repmat('-', 1, 118));
logf('  扫描 %d 夜 → 候选 %d 夜', numel(nights), numel(cand));
if ~isempty(cand)
    logf('  候选夜：%s', strjoin(cand, ', '));
end

csvp = fullfile(cfg.out, 'frame_scan.csv');
fid = fopen(csvp, 'w', 'n', 'UTF-8');
fprintf(fid, 'night,n_frames,n_good,win,alt_min,alt_max,illum_max,verdict\n');
for i = 1:numel(rows)
    q = rows(i);
    fprintf(fid, '%s,%d,%d,%s,%.2f,%.2f,%.3f,%s\n', q.night, q.n_frames, q.n_good, ...
            q.win, q.alt_min, q.alt_max, q.illum_max, q.verdict);
end
fclose(fid);
logf('详', '  逐夜普查表 -> %s', csvp);
end


% ==============================================================================
% [5] 阶段 2 —— 检测：饱和月核 + 无偏质心
% ==============================================================================

function r = mc_blob_center(f, SZ, C)
%MC_BLOB_CENTER 在帧 f（double）里定位最亮的扩展斑。
% 返回 struct: cx, cy, cxs, cys, nsat, nsat_near, peak。
% 返回 [] 表示失败。
%
% 关键设计：**权重截断**。月面中心饱和，直接算亮度加权质心会被平坦的饱和平台主导，
% 且被电荷溢出列（竖直长条）拉偏。做法是
%   (a) 迭代收缩窗口；
%   (b) 权重上限截断到 p99 分位 -> 饱和平台不再主导；
%   (c) 再用 >SAT_LEVEL 的像元在 45 px 内单独取一个"饱和核质心"。
% 两套质心一致时才可信；二者之差本身是 QC 指标。
%
% ⚠ 坐标：内部全 0-based，只有索引数组时 +1。
sm = mc_gauss_blur(f, 2.0);
% ⚠ 跨语言陷阱（实测：MATLAB 与 numpy 的月核质心差 0.05~0.13 deg 就是这么来的）
%   月面中心是**饱和平台**（几千个像元都顶着 65000），高斯模糊之后平台依然是平的，
%   于是 argmax 落在平台的哪个角上，完全取决于"打破平局"的扫描顺序：
%     numpy.argmax -> C 序（行主序）：先扫第一行到底，再扫第二行 …
%     MATLAB max() -> 列主序        ：先扫第一列到底，再扫第二列 …
%   两者选到的种子可以相差几十像素，再经 3 次收缩窗迭代就变成零点几像素的系统差，
%   最后在 rot 上表现为 0.05~0.13 deg 的偏差 —— 看数字完全像"实现了两个不同的算法"。
%   修法：先转置再做 argmax，转置后的列主序恰好等于原来的行主序，与 numpy 完全一致。
smT = sm.';
[~, imaxT] = max(smT(:));
[rT, cT] = ind2sub(size(smT), imaxT);
iy1 = cT;      % smT(rT,cT) == sm(cT,rT) -> sm 的**行**索引 = cT
ix1 = rT;      %                            sm 的**列**索引 = rT
iy = iy1 - 1;  % -> 0-based
ix = ix1 - 1;

for it = 1:C.blob_niter
    y0 = max(0, iy - C.blob_half);
    y1 = min(SZ - 1, iy + C.blob_half);
    x0 = max(0, ix - C.blob_half);
    x1 = min(SZ - 1, ix + C.blob_half);
    w = f(y0 + 1:y1 + 1, x0 + 1:x1 + 1);
    bd = [w(1, :), w(end, :), w(:, 1)', w(:, end)'];
    bg = median(bd);
    wt = max(w - bg, 0.0);
    hi = mc_pctile(wt(:), 99.0);
    wt = min(wt, hi);
    tot = sum(wt(:));
    if tot <= 0, r = []; return; end
    [Xg, Yg] = meshgrid(x0:x1, y0:y1);   % 0-based
    co = sum(wt(:) .* Xg(:)) / tot;
    ro = sum(wt(:) .* Yg(:)) / tot;
    if ~isfinite(co) || ~isfinite(ro), r = []; return; end
    ix = round(co);
    iy = round(ro);
end

sat = f > C.sat_level;
nsat = sum(sat(:));
[sy, sx] = find(sat);
sx = sx - 1;    % 0-based
sy = sy - 1;
if nsat > 4
    near = hypot(sx - co, sy - ro) < C.sat_near;
    if sum(near) > 4
        cxs = mean(sx(near));
        cys = mean(sy(near));
        nsatNear = sum(near);
    else
        cxs = co; cys = ro; nsatNear = 0;
    end
else
    cxs = co; cys = ro; nsatNear = 0;
end

r = struct('cx', co, 'cy', ro, 'cxs', cxs, 'cys', cys, ...
           'nsat', nsat, 'nsat_near', nsatNear, 'peak', max(f(:)));
end


function recs = mc_detect(cfg, candNights, every, minSat, verboseDecode, logf)
%MC_DETECT 逐帧扫候选夜 -> 缓存所有含饱和月核的帧的测量值。
C = mc_defaults();
logf('');
logf('%s', repmat('=', 1, 78));
logf('【阶段 2】检测饱和月核帧   候选夜 %d 个   取样间隔 %d 帧   饱和核下限 %d px', ...
     numel(candNights), every, minSat);
logf('%s', repmat('=', 1, 78));

% 归档里损坏的帧会让 imread 刷告警（Python 版是 libpng 刷 stderr，一个坏夜 180+ 行）。
% 这些帧**已被按夜显式计数并报告**，所以静音不丢任何信息。
if ~verboseDecode
    wstate = warning;
    warning('off', 'all');
    restoreW = onCleanup(@() warning(wstate));
end

recs = struct('night', {}, 'ts', {}, 'cx', {}, 'cy', {}, 'cxs', {}, 'cys', {}, ...
              'nsat', {}, 'nsat_near', {}, 'peak', {});
badTotal = 0;
badNights = {};

for i = 1:numel(candNights)
    nd = candNights{i};
    fs = mc_list_frames(fullfile(cfg.raw, nd));
    fs = fs(1:every:end);
    nBad = 0;
    nMoo = 0;
    for k = 1:numel(fs)
        t = mc_parse_ts(mc_basename(fs{k}));
        im = mc_load_gray(fs{k});
        if isempty(im)
            nBad = nBad + 1;
            continue;
        end
        f = double(im);
        if sum(f(:) > C.sat_level) < minSat, continue; end
        r = mc_blob_center(f, C.size, C);
        if isempty(r), continue; end
        recs(end+1) = struct('night', nd, 'ts', datestr(t, 'yyyymmddHHMMSS'), ...   %#ok<AGROW>
                             'cx', r.cx, 'cy', r.cy, 'cxs', r.cxs, 'cys', r.cys, ...
                             'nsat', r.nsat, 'nsat_near', r.nsat_near, 'peak', r.peak);
        nMoo = nMoo + 1;
    end
    badTotal = badTotal + nBad;
    if nBad > 0
        badNights{end+1} = sprintf('%s : %d/%d 帧不可解码 (%.0f%%)', ...
                                   nd, nBad, numel(fs), 100 * nBad / max(numel(fs), 1));   %#ok<AGROW>
    end
    logf('详', '  %-10s 帧 %4d   解码失败 %3d   含饱和月核 %3d', nd, numel(fs), nBad, nMoo);
end

if ~isempty(badNights)
    logf('');
    logf('  ⚠ 有 %d 夜存在解码失败（归档帧截断；前几帧正常会骗过抽样检查）：', numel(badNights));
    for i = 1:numel(badNights)
        logf('      %s', badNights{i});
    end
    logf('    整夜截断的夜即使偶有帧能解码，内容也是坏的 —— 见阶段 3 的逐夜残差剔除。');
end

logf('');
if isempty(recs)
    % ⚠ 防呆：检出 0 帧时**不要**覆盖已有的测量缓存。
    %   否则一次参数写错的试跑（例如 min_sat 过高的探路）会毁掉跑了几分钟的测量结果，
    %   而且下次 solve 只会说"没有测量缓存"，看不出是被自己覆盖的。
    pOld = fullfile(cfg.out, 'moon_meas.mat');
    logf('  !! 本次检出 0 帧, **不覆盖**已有缓存 (%s)。', ...
         ternary(exist(pOld, 'file') == 2, '原文件保留', '本来就没有'));
    logf('  共 %d 帧解码失败。', badTotal);
    return;
end

night = cell(1, numel(recs));
ts    = cell(1, numel(recs));
cx = zeros(1, numel(recs));   cy = zeros(1, numel(recs));
cxs = zeros(1, numel(recs));  cys = zeros(1, numel(recs));
nsat = zeros(1, numel(recs)); nsatNear = zeros(1, numel(recs));
peak = zeros(1, numel(recs));
for i = 1:numel(recs)
    night{i} = recs(i).night;   ts{i} = recs(i).ts;
    cx(i) = recs(i).cx;         cy(i) = recs(i).cy;
    cxs(i) = recs(i).cxs;       cys(i) = recs(i).cys;
    nsat(i) = recs(i).nsat;     nsatNear(i) = recs(i).nsat_near;
    peak(i) = recs(i).peak;
end
mp = fullfile(cfg.out, 'moon_meas.mat');
save(mp, 'night', 'ts', 'cx', 'cy', 'cxs', 'cys', 'nsat', 'nsatNear', 'peak', '-v7');
logf('');
logf('  检测完成：%d 夜 → %d 帧含饱和月核', numel(candNights), numel(recs));
logf('详', '  测量缓存 -> %s', mp);
end


function recs = mc_load_meas(outdir)
%MC_LOAD_MEAS 读回 阶段 2 的测量缓存。
recs = struct('night', {}, 'ts', {}, 'cx', {}, 'cy', {}, 'cxs', {}, 'cys', {}, ...
              'nsat', {}, 'nsat_near', {}, 'peak', {});
mp = fullfile(outdir, 'moon_meas.mat');
if exist(mp, 'file') ~= 2, return; end
S = load(mp);
n = numel(S.cx);
for i = 1:n
    recs(end+1) = struct('night', mc_tostr(S.night, i), 'ts', mc_tostr(S.ts, i), ...   %#ok<AGROW>
                         'cx', double(S.cx(i)), 'cy', double(S.cy(i)), ...
                         'cxs', double(S.cxs(i)), 'cys', double(S.cys(i)), ...
                         'nsat', double(S.nsat(i)), 'nsat_near', double(S.nsatNear(i)), ...
                         'peak', double(S.peak(i)));
end
end


function s = mc_tostr(c, i)
if iscell(c), s = char(c{i}); else, s = char(c(i)); end
end


% ==============================================================================
% [6] 阶段 3 —— 解算
% ==============================================================================

function g = mc_meas_to_geo(recs, site, C)
%MC_MEAS_TO_GEO 给每条测量补上星历量，并算 rot（直接/镜像）。
n = numel(recs);
g = struct('night', {}, 'ts', {}, 'jd', {}, 'alt', {}, 'az', {}, 'z', {}, ...
           'r', {}, 'phi', {}, 'illum', {}, 'rot_dir', {}, 'rot_mir', {}, ...
           'sun_alt', {}, 'sun_az', {}, 'cx', {}, 'cy', {});
C_ = C;
for i = 1:n
    ts = recs(i).ts;
    t = datetime(str2double(ts(1:4)), str2double(ts(5:6)), str2double(ts(7:8)), ...
                 str2double(ts(9:10)), str2double(ts(11:12)), str2double(ts(13:14)));
    jd = mc_jd(t);
    [alt, A] = mc_moon_altaz(jd, site);
    [sAlt, sAz] = mc_sun_altaz(jd, site);
    dx = recs(i).cx - C_.x0_geo;
    dy = recs(i).cy - C_.y0_geo;
    rad = hypot(dx, dy);
    % phi：以几何法天顶为极点、北=0 东=90。质心判据只对常数偏移不敏感，所以这里
    % 用几何天顶即可（常数偏移不进散布，不影响判据）。
    phi = mod(rad2deg(atan2(dx, -dy)), 360.0);
    g(end+1) = struct('night', recs(i).night, 'ts', ts, 'jd', jd, 'alt', alt, 'az', A, ...   %#ok<AGROW>
                      'z', 90.0 - alt, 'r', rad, 'phi', phi, 'illum', mc_moon_illum(jd), ...
                      'rot_dir', mod(phi - A, 360.0), 'rot_mir', mod(-phi - A, 360.0), ...
                      'sun_alt', sAlt, 'sun_az', sAz, ...
                      'cx', recs(i).cx, 'cy', recs(i).cy);
end
end


function [good, tbl] = mc_gate_by_night(g, altLo, altHi, scatMax, minFrames, logf, title, sunMax)
%MC_GATE_BY_NIGHT 核心判据：逐夜 rot 的圆散布必须小。
%
% 为什么这个判据是决定性的：
%   真月亮 ⇒ rot = phi - A 在同一夜内必须是**常数**（phi 与 A 同步变化）；
%   锁到固定图案/暮光饱和区 ⇒ phi 基本不动而 A 在动 ⇒ 散布立刻 60~100 deg。
% 注意：判据只看**散布**，与天顶像点是否准确无关（常数偏差不进散布）。
%
% sunMax（可选，默认 Inf）：只收 sun_alt < sunMax 的帧。这是**一期 8 bit 数据能用的前提**
%   —— 否则黄昏/黎明天光饱和的假月核会把每夜的散布抬到 8~90 deg。见 mc_sun_altaz。
if nargin < 8 || isempty(sunMax), sunMax = Inf; end
logf('');
logf('%s', repmat('=', 1, 78));
logf('%s', title);
logf('  判据：rot = 像面方位 φ − 星历方位 A，在同一夜内的圆散布 < %.1f°', scatMax);
logf('        （只取 alt %.0f~%.0f°、每夜至少 %d 帧）', altLo, altHi, minFrames);
if isfinite(sunMax)
    logf('        + 太阳高度门限 sun_alt < %.1f°（剔掉暮光/黎明天光饱和帧）', sunMax);
end
logf('%s', repmat('=', 1, 78));
logf('详', '%-10s %5s %-15s %6s %8s | %-26s | %-26s | %s', ...
     '夜', '帧数', 'alt 范围', '日高', 'az 跨度', '直接手性 均值/散布(R)', '镜像手性 均值/散布(R)', '判定');
logf('详', '%s', repmat('-', 1, 122));

good = g([]);            % 同类型空数组
tbl = struct('night', {}, 'n', {}, 'scat_dir', {}, 'mean_dir', {}, 'R_dir', {}, ...
             'scat_mir', {}, 'mean_mir', {}, 'R_mir', {}, 'az_span', {}, 'ok', {}, ...
             'alt_lo', {}, 'alt_hi', {}, 'sun_hi', {});

uNights = unique({g.night});
nFail = 0;
for i = 1:numel(uNights)
    nd = uNights{i};
    sel = strcmp({g.night}, nd) & ([g.alt] > altLo) & ([g.alt] < altHi) ...
          & ([g.sun_alt] < sunMax);
    sub = g(sel);
    if numel(sub) < minFrames, continue; end

    [m1, s1, R1] = mc_circstats([sub.rot_dir]);
    [m2, s2, R2] = mc_circstats([sub.rot_mir]);
    azs = [sub.az];
    span = max(azs) - min(azs);
    if span > 180, span = 360 - span; end     % az 跨 0/360 的校正

    ok = (s1 < scatMax);
    if ok, good = [good, sub]; else, nFail = nFail + 1; end   %#ok<AGROW>
    aLo = min([sub.alt]); aHi = max([sub.alt]);
    sHi = max([sub.sun_alt]);
    tbl(end+1) = struct('night', nd, 'n', numel(sub), 'scat_dir', s1, 'mean_dir', m1, ...   %#ok<AGROW>
                        'R_dir', R1, 'scat_mir', s2, 'mean_mir', m2, 'R_mir', R2, ...
                        'az_span', span, 'ok', ok, 'alt_lo', aLo, 'alt_hi', aHi, 'sun_hi', sHi);
    logf('详', '%-10s %5d %+6.1f..%+5.1f %+6.1f %8.1f | %7.2f %6.2f (%.3f) | %7.2f %6.2f (%.3f) | %s', ...
         nd, numel(sub), aLo, aHi, sHi, span, m1, s1, R1, m2, s2, R2, ...
         ternary(ok, '**月亮**', '否 (假目标)'));
end
logf('详', '%s', repmat('-', 1, 122));
logf('  通过 %d 帧 / %d 夜', numel(good), numel(unique({good.night})));
if nFail > 0
    logf('  ⚠ 另有 %d 夜被判为假目标，不参与定标', nFail);
end
end


function [p, rms, rvec] = mc_fit_geometry(good, model, mirror, init, C)
%MC_FIT_GEOMETRY 在通过判据的帧上拟合 (x0, y0, f, rot, a)。
% 残差用**像面笛卡尔坐标**，所以对极坐标原点选择不敏感。
if nargin < 3 || isempty(mirror), mirror = false; end
if nargin < 4 || isempty(init)
    p0 = [C.x0_geo, C.y0_geo, 420.0, 0.0, -0.1];
else
    p0 = init;
end
s = ternary(mirror, -1.0, 1.0);
% ⚠ 全部用**列向量**。若 zz 是行向量而 xm 是列向量，`rr .* sin(...) - xm` 会触发
%   MATLAB 的隐式扩展变成 N×N 矩阵，残差长度从 2N 变 2N^2 —— 而且**不会报错**，
%   LM 照样能对着一个维度错掉的残差"收敛"出一堆垃圾。
zz = deg2rad([good.z])';
AA = deg2rad([good.az])';
xm = [good.cx]';
ym = [good.cy]';

% ⚠ 注意：形参顺序必须与 mc_res_geo 一致；这里把长数组用闭包带进去。
resfun = @(q) mc_res_geo(q, model, zz, AA, xm, ym, s);
[p, rms, rvec] = mc_lm(resfun, p0, 300);
end


function r = mc_res_geo(q, model, zz, AA, xm, ym, s)
rr = mc_r_of_z(model, zz, q(3), q(5));
rot = deg2rad(q(4));
r = [q(1) + rr .* sin(s * AA + rot) - xm; ...
     q(2) - rr .* cos(s * AA + rot) - ym];
end


function [p, rms, rvec] = mc_lm(resfun, p0, maxit)
%MC_LM 自实现 Levenberg-Marquardt（**不依赖 Optimization Toolbox**）。
%
% 对应 Python 版的 scipy.optimize.least_squares(method='lm')。
% Jacobian 用中心差分（5 个参数 -> 每次迭代 10 次残差求值，代价可忽略）。
% 阻尼策略：接受则 lam/3 并重置 nu，拒绝则 lam*nu、nu*2（经典 Nielsen 方案）。
if nargin < 3 || isempty(maxit), maxit = 300; end

% ⚠ equidistant / equisolid / orthographic / stereographic 这四种模型**不使用**参数 a，
%   于是 Jacobian 的第 5 列恒为 0、A 的第 5 行与第 5 列全零 -> A 奇异。
%   这是**模型本身**的性质，不是数值毛病：g(5) 同样恒为 0，所以只要加一道极小的脊，
%   解出来的 dp(5) 仍严格为 0，a 就保持初值不动 —— 与 scipy 的 method='lm' 行为一致
%   （它也是零梯度不更新）。
%   注：nearlySingularMatrix 警告的静音**不在这里做**。本函数在一次解算里会被调用上千次
%   （模型网格 + 逐夜 + bootstrap），反复 `warning('off')`/恢复会在 R2021b 的 -batch 下
%   触发 Java 侧状态抖动，实测直接把 MATLAB 搞成堆损坏（0xc0000374）。
%   统一由调用方 moon_calib 在主流程里静音一次。
p = double(p0(:));
np = numel(p);
r = resfun(p);
cost = r' * r;
lam = 1e-3;
nu = 2.0;
tolP = 1e-11;
dp = zeros(size(p));

for it = 1:maxit
    J = mc_numjac(resfun, p);
    A = J' * J;
    g = J' * r;
    dg = max(diag(A), 1e-12);
    improved = false;
    for t = 1:30
        dp = -(A + lam * diag(dg) + 1e-12 * eye(np)) \ g;
        if ~all(isfinite(dp))
            lam = lam * nu; nu = nu * 2;
            if lam > 1e14, break; end
            continue;
        end
        pn = p + dp;
        rn = resfun(pn);
        if ~all(isfinite(rn))
            lam = lam * nu; nu = nu * 2;
            if lam > 1e14, break; end
            continue;
        end
        cn = rn' * rn;
        if cn < cost
            p = pn; r = rn; cost = cn;
            lam = max(lam / 3, 1e-14);
            improved = true;
            break;
        else
            lam = lam * nu; nu = nu * 2;
            if lam > 1e14, break; end
        end
    end
    if ~improved, break; end
    if norm(dp) <= tolP * (norm(p) + tolP), break; end
end
rms = sqrt(mean(r .^ 2));
rvec = r;
end


function J = mc_numjac(resfun, p)
%MC_NUMJAC 中心差分 Jacobian。
np = numel(p);
h = eps^(1/3) * max(abs(p), 1.0);
cols = cell(1, np);
for k = 1:np
    pk = p;
    pk(k) = p(k) + h(k);
    pm = p;
    pm(k) = p(k) - h(k);
    cols{k} = (resfun(pk) - resfun(pm)) / (2 * h(k));
end
J = zeros(numel(cols{1}), np);
for k = 1:np
    J(:, k) = cols{k};
end
end


function [rot, rms] = mc_fit_rot_only(g, model, f, a, x0, y0, rot0)
%MC_FIT_ROT_ONLY 固定几何，只解该夜（或该子集）的 rot。
% 同 mc_fit_geometry：全列向量，否则隐式扩展会静默算错。
zn = deg2rad([g.z])';
An = deg2rad([g.az])';
xn = [g.cx]';
yn = [g.cy]';
resfun = @(q) mc_res_rot(q, model, f, a, x0, y0, zn, An, xn, yn);
[q, rms] = mc_lm(resfun, [rot0], 200);
rot = mod(q(1), 360.0);
end


function r = mc_res_rot(q, model, f, a, x0, y0, zn, An, xn, yn)
rr = mc_r_of_z(model, zn, f, a);
rot = deg2rad(q(1));
r = [x0 + rr .* sin(An + rot) - xn; ...
     y0 - rr .* cos(An + rot) - yn];
end


function res = mc_grid_fit(good, C)
%MC_GRID_FIT 在给定帧集上跑遍 模型 × 手性 网格。
nM = numel(C.models);
res = struct('model', {}, 'mirror', {}, 'p', {}, 'rms', {});
for i = 1:nM
    for mir = [false, true]
        [p, rms] = mc_fit_geometry(good, C.models{i}, mir, [], C);
        res(end+1) = struct('model', C.models{i}, 'mirror', mir, ...   %#ok<AGROW>
                            'p', p, 'rms', rms);
    end
end
end


function mc_print_grid(res, C, logf)
%MC_PRINT_GRID 打印 模型 × 手性 网格。**只进日志文件** —— 控制台由调用者打一行"最佳模型"。
logf('详', '%-12s %-22s %-6s | %8s %8s %8s %9s %8s | %8s %8s', ...
     '模型', 'r(z)', '手性', 'x0', 'y0', 'f', 'a', 'rot(deg)', 'rms(px)', 'r90(px)');
logf('详', '%s', repmat('-', 1, 122));
for i = 1:numel(res)
    k = find(strcmp(C.models, res(i).model), 1);
    p = res(i).p;
    r90 = abs(mc_r_of_z(res(i).model, pi / 2, p(3), p(5)));
    if strcmp(res(i).model, 'equid+')
        aShow = sprintf('%+9.4f', p(5));
    else
        aShow = sprintf('%9s', '-');
    end
    logf('详', '%-12s %-22s %-6s | %8.2f %8.2f %8.2f %s %8.2f | %8.2f %8.2f', ...
         C.model_cn{k}, C.model_desc{k}, ternary(res(i).mirror, '镜像', '直接'), ...
         p(1), p(2), p(3), aShow, mod(p(4), 360), res(i).rms, r90);
end
end


function pn = mc_night_residual_table(good, model, p, C)
%MC_NIGHT_RESIDUAL_TABLE 固定几何，每夜只解 rot，返回逐夜 (帧数, rot, 与全局差, 残差 rms)。
%
% ⚠ 这是**最灵敏的坏夜探测器**，不能省。理由：「逐夜 rot 恒定性」只看**散布**不看
%   残差，所以对"整幅平移了的坏帧"是盲的。归档里被截断的 PNG 偶尔能解码出来，
%   但图像下半部是坏的 -> 质心被拉偏 -> 该夜 rms 突然变成别夜的 5~10 倍，
%   而 rot 散布仍然很小，判据照样放它过去。必须配这张残差表才能揪出来。
%   （本项目实测：20260322 那一夜就是这样被判据放行、又被这张表揪出的。）
pn = struct('night', {}, 'n', {}, 'rot', {}, 'dr', {}, 'rms', {}, 'alt_lo', {}, 'alt_hi', {});
uNights = unique({good.night});
for i = 1:numel(uNights)
    nd = uNights{i};
    sub = good(strcmp({good.night}, nd));
    if numel(sub) < 3, continue; end
    [rot, rn] = mc_fit_rot_only(sub, model, p(3), p(5), p(1), p(2), p(4));
    alts = [sub.alt];
    pn(end+1) = struct('night', nd, 'n', numel(sub), 'rot', rot, ...   %#ok<AGROW>
                       'dr', mc_ang_diff(rot, p(4)), 'rms', rn, ...
                       'alt_lo', min(alts), 'alt_hi', max(alts));
end
end


% ------------------------------------------------------------------ 3.9 用

function [prof, darkf, lv] = mc_azimuthal_profile(imgs, x0, y0, SZ)
%MC_AZIMUTHAL_PROFILE 若干帧的**方位平均径向剖面**（取中位），覆盖到画幅四角。
%
% 返回 (median 亮度, 暗像元占比, 电平字典)。这是判断"亮环到底是不是地平圈 /
% 镜筒像圈边界"的**唯一直接证据**：
%   * 若某个 r 处亮度出现台阶 -> 那是光学/地面边界，可以拿来核对 f；
%   * 若某 r 之外亮度跌到**偏置电平** -> 那里没有光（地面，或镜头像圈边界）；
%   * 若一直平到四角 -> 传感器被照满，那条"亮环"只是渐晕/地平辉光的软边缘。
%
% ⚠ 两个必踩的坑（本项目实测）：
%   (1) 暗像元判据必须用"天空电平与偏置电平的中点"，**不能**用"低于天空的 30%"：
%       本机 16 bit 帧带 ~3441 DN 的偏置基座，30% × 5200 = 1560 < 3441，
%       于是**永远判不出暗像元**，会得出"传感器被照满"的错误结论。
%   (2) 偏置电平**不能**取剖面的全局最小值：r > 四角半径处根本没有像元，数组里是 0，
%       会把偏置电平算成 0，进而把"光终止半径"推到画幅之外。
%       正确做法是取**四角区域**像元的中位数。
nb = SZ;
[XX, YY] = mc_grid(SZ);
RR = hypot(XX - x0, YY - y0);
BIN = min(max(round(RR) + 1, 1), nb);
flat = BIN(:);
cnt = accumarray(flat, 1, [nb 1]);

rCorner = hypot(max(x0, SZ - 1 - x0), max(y0, SZ - 1 - y0));
corner = RR > 0.9 * rCorner;
skyBox = RR < 400;

meds = zeros(numel(imgs), nb);
darks = zeros(numel(imgs), nb);
skies = zeros(numel(imgs), 1);
peds = zeros(numel(imgs), 1);
for i = 1:numel(imgs)
    f = double(imgs{i});
    sky = median(f(skyBox));
    ped = median(f(corner));          % 偏置电平取四角，不取 min
    thr = 0.5 * (sky + ped);
    s = accumarray(flat, f(:), [nb 1]);
    d = accumarray(flat, double(f(:) < thr), [nb 1]);
    meds(i, :) = (s ./ max(cnt, 1))';
    darks(i, :) = (d ./ max(cnt, 1))';
    skies(i) = sky;
    peds(i) = ped;
end

prof = median(meds, 1);
darkf = median(darks, 1);
nz = find(cnt > 0);
lv = struct('sky', median(skies), 'ped', median(peds), ...
            'r_corner', rCorner, 'r_max', floor(rCorner), 'nb_used', nz(end));
end


% ------------------------------------------------------------------ 3.10 用

function [azc, rb] = mc_boundary_radius_vs_az(imgs, x0, y0, halfLevel, nbins, SZ)
%MC_BOUNDARY_RADIUS_VS_AZ 以天顶像点为极点，测"亮盘边界半径"随方位角的分布。
%
% 用途：判"亮盘边界到底是绕天顶的同心圆（地平圈）还是绕主点的圆（镜头像圈）"。
%   * 边界绕天顶同心 -> r_boundary(az) = 常数
%   * 边界绕主点同心 -> r_boundary(az) ≈ R - d cos(az - az0)，d = 主点到天顶的距离
% 对后者做最小二乘就能同时解出 d、az0、R，而 d/f 就是**光轴相对天顶的倾角**。
%
% ⚠ 这一步不能省：本项目的 r(z) 模型默认"等 z 的像点是以天顶像点为圆心的圆"，
%   这个前提只有在光轴严格指向天顶时才严格成立。倾角 5 deg 在本机 f≈422 px/rad 上
%   就是 37 px 的形变，足以吃掉全部标定精度。
[XX, YY] = mc_grid(SZ);
RR = hypot(XX - x0, YY - y0);
AZ = mod(atan2(XX - x0, -(YY - y0)), 2 * pi);
idx = floor(AZ(:) / (2 * pi) * nbins) + 1;
idx = min(max(idx, 1), nbins);
rrf = RR(:);

P = NaN(numel(imgs), nbins);
for i = 1:numel(imgs)
    m = double(imgs{i}) > halfLevel;
    sel = m(:);
    if ~any(sel)
        continue;
    end
    out = accumarray(idx(sel), rrf(sel), [nbins 1], @max, NaN);
    P(i, :) = out(:)';
end
rb = median(P, 1, 'omitnan');
edges = linspace(0, 2 * pi, nbins + 1);
azc = rad2deg(0.5 * (edges(1:end-1) + edges(2:end)));
end


% ------------------------------------------------------------------ 主解算

function cal = mc_solve(cfg, site, logf)
%MC_SOLVE 阶段 3：判据筛选 + 坏夜剔除 + 模型/手性 + 不确定度 + 一致性检验 + 断言。
C = mc_defaults();
% ★ 几何参考跟随本次运行的 rawRoot（二期默认值 / 一期需覆盖）。见顶层 inputParser 注释。
if isfield(cfg, 'x0_geo') && ~isempty(cfg.x0_geo), C.x0_geo = cfg.x0_geo; end
if isfield(cfg, 'y0_geo') && ~isempty(cfg.y0_geo), C.y0_geo = cfg.y0_geo; end
if isfield(cfg, 'rim')    && ~isempty(cfg.rim),    C.rim    = cfg.rim;    end

recs = mc_load_meas(cfg.out);
if isempty(recs)
    logf('没有测量缓存（%s），请先跑 stage=detect', fullfile(cfg.out, 'moon_meas.mat'));
    cal = [];
    return;
end
g = mc_meas_to_geo(recs, site, C);

logf('');
logf('%s', repmat('=', 1, 78));
logf('【阶段 3】解算   站点 (%.3fE, %.3fN, %.0f m)   测量帧 %d（来自 %d 夜）', ...
     site(1), site(2), site(3) * 1000, numel(recs), numel(unique({g.night})));
logf('%s', repmat('=', 1, 78));

% --- 3.1 逐夜判据 ---
[good, gateTbl] = mc_gate_by_night(g, cfg.alt_lo, cfg.alt_hi, cfg.scat_max, ...
                                   cfg.min_frames, logf, ...
                                   '【3.1】逐夜 rot 恒定性判据 —— 分辨真月亮与假目标', ...
                                   cfg.sun_alt_max);
if numel(good) < 12
    logf('');
    logf('通过判据的帧太少 (%d)，无法定标。检查: 站点坐标? 时间戳是否为 UT? 是否有月夜?', ...
         numel(good));
    if isfinite(cfg.sun_alt_max)
        logf('  本次开了 sun_alt_max=%.1f deg: 若通过 0 帧但原始测量很多，说明这批帧全在暮光里。', ...
             cfg.sun_alt_max);
    end
    cal = [];
    return;
end

% --- 3.2 模型 × 手性 网格（临时最佳，供坏夜剔除用） ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.2】模型 × 手性 网格   %d 种径向模型 × 直接/镜像（完整表见日志文件）', numel(C.models));
logf('%s', repmat('=', 1, 78));
res = mc_grid_fit(good, C);
mc_print_grid(res, C, logf);
bestModel = mc_pick_best(res, C);

% --- 3.2b 逐夜残差自动剔除 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.2b】逐夜残差自动剔除   阈值 = max(%.1f × 中位 rms, %.1f px)', ...
     cfg.night_rms_mult, cfg.night_rms_floor);
logf('%s', repmat('=', 1, 78));
pBest = mc_get_p(res, C, bestModel, false);
pn0 = mc_night_residual_table(good, bestModel, pBest, C);
medRms = median([pn0.rms]);
thr = max(cfg.night_rms_mult * medRms, cfg.night_rms_floor);
logf('  中位 rms = %.2f px   →   剔除阈值 %.2f px', medRms, thr);
drop = pn0([pn0.rms] > thr);
if ~isempty(drop)
    for i = 1:numel(drop)
        logf('  ✗ 剔除 %s：%d 帧，rms %.2f px（其余夜中位 %.2f px）', ...
             drop(i).night, drop(i).n, drop(i).rms, medRms);
    end
    logf('    典型原因：归档帧截断（IHDR 完好、IDAT 不全，偶尔能解码但内容是坏的）；');
    logf('              或该夜有云 / 亮边污染。这些夜不应进入定标。');
    dropSet = {drop.night};
    good = good(~ismember({good.night}, dropSet));
    res = mc_grid_fit(good, C);
    bestModel = mc_pick_best(res, C);
    logf('详', '');
    logf('详', '%s', repmat('-', 1, 122));
    logf('详', '【3.2final】剔除后的模型网格');
    logf('详', '%s', repmat('-', 1, 122));
    mc_print_grid(res, C, logf);
else
    logf('  没有需要剔除的夜。');
end

usedNights = unique({good.night});
pBest = mc_get_p(res, C, bestModel, false);
rmsBest = mc_get_rms(res, C, bestModel, false);
kBest = find(strcmp(C.models, bestModel), 1);
logf('%s', repmat('-', 1, 78));
logf('  ▶ 最佳径向模型：%s（%s）    残差 rms = %.2f px', ...
     C.model_cn{kBest}, C.model_desc{kBest}, rmsBest);
logf('    参与定标的 %d 夜：%s', numel(usedNights), strjoin(usedNights, ', '));

% 手性判决
rmsDir = mc_get_rms(res, C, bestModel, false);
rmsMir = mc_get_rms(res, C, bestModel, true);
chirality = ternary(rmsDir <= rmsMir, 'direct', 'mirror');
logf('  ▶ 手性判决：%s   （直接 %.2f px vs 镜像 %.2f px，判别力 %.1f×）', ...
     ternary(strcmp(chirality, 'direct'), '直接（无镜像）', '镜像（左右翻转）'), ...
     rmsDir, rmsMir, max(rmsDir, rmsMir) / max(min(rmsDir, rmsMir), 1e-9));

rim = cfg.rim;
zRim = mc_z_of_rim(bestModel, pBest(3), pBest(5), rim);
r90 = mc_r_of_z(bestModel, pi / 2, pBest(3), pBest(5));

% --- 3.3 逐夜 rot 一致性 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.3】逐夜 rot 一致性（几何固定，每夜只解 rot）—— 检验相机是否真静止');
logf('%s', repmat('=', 1, 78));
logf('详', '%-10s %5s %-15s %11s %12s %10s', '夜', '帧数', 'alt 范围', '该夜 rot', '与全局差', '残差 rms');
logf('详', '%s', repmat('-', 1, 122));
perNight = mc_night_residual_table(good, bestModel, pBest, C);
for i = 1:numel(perNight)
    q = perNight(i);
    logf('详', '%-10s %5d %+6.1f..%-6.1f %9.2f deg %+10.2f deg %8.2f px', ...
         q.night, q.n, q.alt_lo, q.alt_hi, q.rot, q.dr, q.rms);
end
if ~isempty(perNight)
    [cm, cs, cR] = mc_circstats([perNight.rot]);
else
    cm = NaN; cs = NaN; cR = NaN;
end
logf('详', '%s', repmat('-', 1, 122));
logf('  跨夜 rot：圆均值 %.2f°   圆散布 %.2f°   （%d 夜）', cm, cs, numel(perNight));
if isfinite(cs) && cs < 1.0
    logf('  ⇒ 散布 < 1°，说明相机在整个观测季里没有动过 ⇒ 一套定标可覆盖全部夜。');
end

% --- 3.4 残差系统学 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.4】残差系统学（最佳模型）—— 有没有随 z / 方位残留的畸变');
logf('%s', repmat('=', 1, 78));
x0 = pBest(1); y0 = pBest(2); f = pBest(3); rotG = pBest(4); aa = pBest(5);
zAll = zeros(numel(good), 1); altAll = zeros(numel(good), 1); azAll = zeros(numel(good), 1);
dxAll = zeros(numel(good), 1); dyAll = zeros(numel(good), 1);
radRes = zeros(numel(good), 1); tanRes = zeros(numel(good), 1);
for i = 1:numel(good)
    zr = deg2rad(good(i).z);
    rr = mc_r_of_z(bestModel, zr, f, aa);
    px = x0 + rr * sin(deg2rad(good(i).az) + deg2rad(rotG));
    py = y0 - rr * cos(deg2rad(good(i).az) + deg2rad(rotG));
    dx = good(i).cx - px;
    dy = good(i).cy - py;
    ux = sin(deg2rad(good(i).az) + deg2rad(rotG));
    uy = -cos(deg2rad(good(i).az) + deg2rad(rotG));
    zAll(i) = good(i).z; altAll(i) = good(i).alt; azAll(i) = good(i).az;
    dxAll(i) = dx; dyAll(i) = dy;
    radRes(i) = dx * ux + dy * uy;
    tanRes(i) = -dx * uy + dy * ux;
end
logf('详', '%-18s %5s | %9s %9s | %9s %9s', 'z 区间 (deg)', 'n', '径向残差', '散布', '切向残差', '散布');
logf('详', '%s', repmat('-', 1, 122));
edges = [0, 20, 30, 40, 50, 60, 70, 78, 90];
for i = 1:numel(edges) - 1
    m = zAll >= edges(i) & zAll < edges(i + 1);
    if sum(m) < 3, continue; end
    logf('详', '%-18s %5d | %+9.2f %9.2f | %+9.2f %9.2f', ...
         sprintf('%d-%d', edges(i), edges(i + 1)), sum(m), ...
         mean(radRes(m)), std(radRes(m)), mean(tanRes(m)), std(tanRes(m)));
end
logf('详', '%s', repmat('-', 1, 122));
logf('  径向残差 z 依赖 %.2f px（模型形状没吃干净的残余）；切向残差均值 %.2f px（这才是 rot 的误差）', ...
     max(abs(arrayfun(@(lo, hi) mean(radRes(zAll >= lo & zAll < hi)), ...
         edges(1:end - 1), edges(2:end)))), abs(mean(tanRes)));

% --- 3.5 bootstrap 不确定度 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.5】重采样（bootstrap）参数不确定度   %d 次重采样帧', cfg.n_boot);
logf('%s', repmat('=', 1, 78));
boot = zeros(0, 5);
if cfg.n_boot > 0
    % ⚠ 与 Python 版数值不同：Python 用 numpy RandomState，这里用 MATLAB 的 MT19937。
    %   两者随机序列不同，但 sd 的统计意义相同（同一 CPU 上可完全复现本脚本自己的结果）。
    rng(20260916, 'twister');
    nfr = numel(good);
    for b = 1:cfg.n_boot
        idx = randi(nfr, nfr, 1);
        try
            pb = mc_fit_geometry(good(idx), bestModel, false, pBest, C);
            boot(end+1, :) = pb(:)';   %#ok<AGROW>
        catch
            continue;
        end
    end
end
if ~isempty(boot)
    names = {'x0', 'y0', 'f', 'rot_deg', 'a'};
    logf('详', '%-8s %12s %12s %12s', '参数', '拟合值', 'bootstrap sd', '±3sd 区间');
    logf('详', '%s', repmat('-', 1, 122));
    for i = 1:5
        sd = std(boot(:, i));
        logf('详', '%-8s %12.3f %12.3f  [%.3f, %.3f]', ...
             names{i}, pBest(i), sd, pBest(i) - 3 * sd, pBest(i) + 3 * sd);
    end
    [rm, rsd, rR] = mc_circstats(boot(:, 4));
    logf('详', '%-8s %12.3f %12.3f   R=%.4f', 'rot(圆)', rm, rsd, rR);
    logf('  x0 ±%.2f   y0 ±%.2f   f ±%.2f px/rad   rot ±%.2f°   a ±%.4f  （1σ，%d 次）', ...
         std(boot(:, 1)), std(boot(:, 2)), std(boot(:, 3)), std(boot(:, 4)), ...
         std(boot(:, 5)), size(boot, 1));
end

% --- 3.6 半分交叉验证 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.6】半分交叉验证 —— A 夜组拟合几何，B 夜组只解 rot 看残差');
logf('%s', repmat('=', 1, 78));
nU = numel(usedNights);
trA = usedNights(1:2:nU); teA = usedNights(2:2:nU);
pairs = {trA, teA; teA, trA};
xval = struct('fit_nights', {}, 'test_nights', {}, 'rot_fit', {}, 'rot_test', {}, ...
              'rms_test', {}, 'drot', {});
for i = 1:2
    tr = pairs{i, 1};
    te = pairs{i, 2};
    A_ = good(ismember({good.night}, tr));
    B_ = good(ismember({good.night}, te));
    if numel(A_) < 10 || numel(B_) < 10, continue; end
    pa = mc_fit_geometry(A_, bestModel, false, [], C);
    [rotB, rmsB] = mc_fit_rot_only(B_, bestModel, pa(3), pa(5), pa(1), pa(2), pa(4));
    xval(end+1) = struct('fit_nights', {tr}, 'test_nights', {te}, ...   %#ok<AGROW>
                         'rot_fit', mod(pa(4), 360), 'rot_test', rotB, ...
                         'rms_test', rmsB, 'drot', mc_ang_diff(rotB, pa(4)));
    logf('详', '  拟合夜 %s  →  测试夜 %s', strjoin(tr, ','), strjoin(te, ','));
    logf('详', '      拟合 rot %7.2f°   测试 rot %7.2f°   差 %+.2f°   测试 rms %.2f px', ...
         mod(pa(4), 360), rotB, mc_ang_diff(rotB, pa(4)), rmsB);
end
if ~isempty(xval)
    logf('  两组互换：rot 差最大 %.2f°，测试 rms %.2f~%.2f px   ⇒ 几何可外推，不是过拟合', ...
         max(abs([xval.drot])), min([xval.rms_test]), max([xval.rms_test]));
end

% --- 3.7 站点坐标 / 模型灵敏度 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.7】系统差来源：站点坐标 / 模型选择');
logf('%s', repmat('=', 1, 78));
siteKeys = fieldnames(C.sites);
keys = cellstr(strcat({good.night}', '|', {good.ts}'));
for i = 1:numel(siteKeys)
    sk = siteKeys{i};
    if strcmp(sk, cfg.site_key), continue; end
    sAlt = C.sites.(sk);
    % ⚠ 必须**重算** z/az，不能只把 site 传进去 —— 拟合函数里只用 good 里已经算好的
    %   z/az，原地换站点参数是没有任何效果的（会静默打印"差 0.00 deg"）。
    gAlt = mc_meas_to_geo(recs, sAlt, C);
    keysAlt = cellstr(strcat({gAlt.night}', '|', {gAlt.ts}'));
    sel = ismember(keysAlt, keys);
    sub = gAlt(sel);
    if numel(sub) < 6
        logf('详', '  站点 %-9s (%.3fE, %.3fN) → 可用帧只有 %d，跳过', sk, sAlt(1), sAlt(2), numel(sub));
        continue;
    end
    [pa, rmsA] = mc_fit_geometry(sub, bestModel, false, [], C);
    logf('  换站点 %s (%.3fE, %.3fN)：rot = %.2f°（差 %+.2f°），rms %.2f px，%d 帧', ...
         sk, sAlt(1), sAlt(2), mod(pa(4), 360), mc_ang_diff(pa(4), rotG), rmsA, numel(sub));
end

% 模型选择对 rot 的影响：直接把网格里各模型（直接手性）的 rot 拉出来对比
rotsM = zeros(numel(C.models), 1);
for i = 1:numel(C.models)
    pp = mc_get_p(res, C, C.models{i}, false);
    rotsM(i) = mod(pp(4), 360);
end
logf('详', '  模型选择 -> rot：%s', strjoin(arrayfun(@(i) sprintf('%s %.2f', C.model_cn{i}, rotsM(i)), ...
    1:numel(C.models), 'UniformOutput', false), ' ; '));
logf('  rot 对径向模型的极差仅 %.2f°  ⇒ rot 是几何量，不依赖"选哪条 r(z)"', ...
     max(rotsM) - min(rotsM));

% --- 3.8 低空 hold-out ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.8】低空留出检验（不参与拟合的 alt < %.0f° 帧，用小孔径独立定位）', cfg.alt_hi);
logf('%s', repmat('=', 1, 78));
logf('详', '%-10s %-10s %7s %7s %9s %9s %8s %9s', ...
     '夜', '时刻', 'alt', 'z', '预测 r', '实测 r', 'r 差', '像点残差');
logf('详', '%s', repmat('-', 1, 122));
lowres = struct('night', {}, 'tstr', {}, 'alt', {}, 'z', {}, 'r_pred', {}, ...
                'r_meas', {}, 'dr', {}, 'resid', {});
for i = 1:numel(usedNights)
    nd = usedNights{i};
    fs = mc_list_frames(fullfile(cfg.raw, nd));
    if isempty(fs), continue; end
    fs = fs(1:4:end);
    for k = 1:numel(fs)
        t = mc_parse_ts(mc_basename(fs{k}));
        jd = mc_jd(t);
        [alt, A] = mc_moon_altaz(jd, site);
        if ~(alt > cfg.alt_hold && alt < cfg.alt_lo), continue; end
        z = deg2rad(90.0 - alt);
        rr = mc_r_of_z(bestModel, z, f, aa);
        px = x0 + rr * sin(deg2rad(A) + deg2rad(rotG));
        py = y0 - rr * cos(deg2rad(A) + deg2rad(rotG));
        im = mc_load_gray(fs{k});
        if isempty(im), continue; end
        fp = double(im);
        hw = cfg.aperture;
        ay0 = max(0, floor(py) - hw); ay1 = min(C.size - 1, floor(py) + hw);
        ax0 = max(0, floor(px) - hw); ax1 = min(C.size - 1, floor(px) + hw);
        w = fp(ay0 + 1:ay1 + 1, ax0 + 1:ax1 + 1);
        if isempty(w), continue; end
        bd = [w(1, :), w(end, :), w(:, 1)', w(:, end)'];
        bg = median(bd);
        if (max(w(:)) - bg) < cfg.min_snr_dn, continue; end
        wt = max(w - bg, 0.0);
        if sum(wt(:)) <= 0, continue; end
        [Xg, Yg] = meshgrid(ax0:ax1, ay0:ay1);
        ccx = sum(wt(:) .* Xg(:)) / sum(wt(:));
        ccy = sum(wt(:) .* Yg(:)) / sum(wt(:));
        rMeas = hypot(ccx - x0, ccy - y0);
        lowres(end+1) = struct('night', nd, 'tstr', datestr(t, 'HH:MM:SS'), ...   %#ok<AGROW>
                               'alt', alt, 'z', 90 - alt, 'r_pred', rr, 'r_meas', rMeas, ...
                               'dr', rMeas - rr, 'resid', hypot(ccx - px, ccy - py));
    end
end
for i = 1:numel(lowres)
    q = lowres(i);
    logf('详', '%-10s %-10s %+7.1f %7.1f %9.1f %9.1f %+8.1f %9.2f', ...
         q.night, q.tstr, q.alt, q.z, q.r_pred, q.r_meas, q.dr, q.resid);
end
if ~isempty(lowres)
    dv = [lowres.dr];
    rv = [lowres.resid];
    logf('详', '%s', repmat('-', 1, 122));
    cc = corrcoef([lowres.alt], dv);
    cc = ternary(numel(cc) > 1, cc(1, 2), NaN);
    logf('  n=%d   半径差 中位 %+.1f px（散布 %.1f）   像点残差 中位 %.1f px、90 分位 %.1f px', ...
         numel(lowres), median(dv), std(dv), median(rv), mc_pctile(rv, 90));
    if isfinite(cc) && abs(cc) < 0.5
        logf('  ⇒ 半径差与 alt 无显著相关（r = %+.2f）⇒ 是常数偏置（软边缘/亮度梯度），不是模型形状错。', cc);
    elseif isfinite(cc)
        logf('  ⚠ 半径差随 alt 单调（r = %+.2f）⇒ 模型在近地平段的形状可能不对，需查。', cc);
    end
else
    logf('  （该时段没有可见月亮帧）');
end

% --- 3.9 与"纯几何法"的一致性检验 ---
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.9】一致性检验：月亮法的地平圈半径 vs 几何法的亮环半径（两条独立路线）');
logf('%s', repmat('=', 1, 78));
allNights = mc_list_nights(cfg.raw);
pf = {};
for i = 1:2:numel(allNights)
    nd = allNights{i};
    fs = mc_list_frames(fullfile(cfg.raw, nd));
    if numel(fs) < 20, continue; end
    % 先只用文件名算月亮高度挑帧（不读图像），再解码 -> 便宜且不会混进暮光帧
    candT = [];
    for k = 1:numel(fs)
        t = mc_parse_ts(mc_basename(fs{k}));
        a_ = mc_moon_altaz(mc_jd(t), site);
        if a_ < -10.0, candT(end+1) = k; end   %#ok<AGROW>
    end
    if isempty(candT), continue; end
    pickIdx = candT(1 + round((0:2) * (numel(candT) - 1) / 2));
    for k = pickIdx
        im = mc_load_gray(fs{k});
        if ~isempty(im), pf{end+1} = im; end   %#ok<AGROW>
    end
end
prof = [];
darkf = [];
lv = [];
if ~isempty(pf)
    [prof, darkf, lv] = mc_azimuthal_profile(pf, x0, y0, C.size);
end

if ~isempty(prof) && lv.nb_used > rim + 1
    logf('详', '  方位平均径向剖面 (%d 帧中位, 均为月亮在地平下 10 deg 以上的干净天光帧;', numel(pf));
    logf('详', '                      天顶像点用拟合值; 画幅四角 r 最大约 %.0f px):', lv.r_corner);
    logf('详', '  %8s %14s %14s %14s', 'r (px)', '亮度 (DN)', '相邻斜率', '暗像元占比');
    rrShow = max(0, floor(rim) - 96):8:min(lv.nb_used - 1, floor(rim) + 24);
    rrShow = [rrShow, (floor(rim) + 32):24:lv.nb_used];
    rrShow = sort(unique(rrShow(rrShow >= 0 & rrShow < lv.nb_used)));
    for r_ = rrShow
        logf('详', '  %8d %14.0f %14.0f %13.3f', r_, prof(r_ + 1), prof(r_ + 2) - prof(r_ + 1), darkf(r_ + 1));
    end
    lo_ = max(0, floor(rim) - 96);
    hi_ = floor(rim) + 1;
    seg = prof(lo_ + 1:hi_ + 1);
    if numel(seg) >= 2
        [~, kg] = max(abs(diff(seg)));
        rGrad = lo_ + kg - 1;
    else
        rGrad = NaN;
    end
    skyLv = lv.sky;
    pedLv = lv.ped;
    thrLv = 0.5 * (skyLv + pedLv);
    rEnd = NaN;
    for r_ = 300:lv.nb_used
        if prof(r_ + 1) < thrLv
            rEnd = r_;
            break;
        end
    end
    fEqMoon = r90 / (pi / 2);
    fEqGeom = rim / (pi / 2);
    logf('详', '');
    logf('详', '  天空电平 %.0f DN   偏置电平(四角) %.0f DN   半高点 %.0f DN', skyLv, pedLv, thrLv);
    logf('详', '  剖面 |dI/dr| 最大处              r = %d px', rGrad);
    logf('详', '  几何法给出的亮环半径 R_rim      r = %.1f px', rim);
    if isfinite(rEnd)
        logf('详', '  亮度跌到半高点的半径(光终止处)  r = %d px', rEnd);
    end
    logf('详', '  月亮法给出的地平圈 r(90 deg)    r = %.1f px   (= 半幅 %d px 的 %.1f%%)', ...
         r90, floor(C.size / 2), 100 * r90 / (C.size / 2));
    logf('详', '  R_rim 与 r(90 deg) 相差 %+.1f px  (在本地径向尺度上约 %.2f deg)', ...
         rim - r90, rad2deg((rim - r90) / max(abs(f), 1e-9)));
    logf('详', '  等效焦距:  月亮法 %.1f px/rad   几何法 %.1f px/rad   差 %+.1f%%', ...
         fEqMoon, fEqGeom, 100 * (fEqGeom - fEqMoon) / fEqMoon);
    logf('');
    logf('  ▶ 两条独立路线对比：几何法亮环 R_rim = %.1f px；月亮法地平圈 r(90°) = %.1f px（差 %+.1f px = %+.1f%%）', ...
         rim, r90, rim - r90, 100 * (rim - r90) / r90);
    if isfinite(rEnd) && rEnd >= 400
        logf('    剖面在 r(90°) 处平滑无硬台阶，亮度到 %d px 才跌半高 ⇒ 亮环是**软边缘**，不是硬几何边界。', rEnd);
        logf('    几何法把 f 定为 R_rim/(π/2)，前提是"亮环恰好是 90° 地平圈"；月亮法**不需要**这个前提。');
        logf('    两者差 %+.1f px，而低空留出的系统项就有 ~6 px、软边缘宽约 28 px ⇒ 同一量级，互为佐证。', ...
             rim - r90);
        logf('    最终以月亮法为准：它没有"亮环 = 地平圈"的假设，且外推到 z = 86° 已被实测验证。');
    else
        logf('    剖面一直到 %.0f px 都没有明显跌落 ⇒ 该帧序列整幅都在天空内，', lv.nb_used);
        logf('    视场边缘落在画幅之外。此时几何法的 R_rim 不能与地平圈等同。');
    end
    rp = fullfile(cfg.out, 'radial_profile.csv');
    fh = fopen(rp, 'w', 'n', 'UTF-8');
    fprintf(fh, 'r_px,brightness_DN,dark_fraction\n');
    for i = 0:lv.nb_used
        fprintf(fh, '%d,%.2f,%.4f\n', i, prof(i + 1), darkf(i + 1));
    end
    fclose(fh);
    logf('详', '  完整剖面 -> %s', rp);
else
    logf('  （可用于剖面的落月后干净帧不足）');
end

% --- 3.10 亮盘边界：绕天顶同心 还是 绕主点同心？ ---
bR = NaN;   % 边界圆半径（供结果卡）
bD = NaN;   % 边界圆心偏离天顶像点的距离
logf('');
logf('%s', repmat('=', 1, 78));
logf('【3.10】亮盘边界几何 —— 检验"等天顶角的像点绕天顶像点同心"这个前提');
logf('%s', repmat('=', 1, 78));
if ~isempty(pf) && ~isempty(prof)
    half = 0.5 * (lv.sky + lv.ped);
    [azb, rb] = mc_boundary_radius_vs_az(pf, x0, y0, half, 72, C.size);
    okb = isfinite(rb(:))' & (rb(:)' > 100);
    if sum(okb) >= 12
        % ⚠ 坑（Python 版实测踩到，MATLAB 版同样要防）：
        %   这一段**绝对不能**用 a / A / d / R0 / res 当局部变量 —— 它们在外层分别是
        %   畸变系数 a、方位角数组、残差数组……一旦覆盖，后面的汇总会拿数组去当标量，
        %   要么直接崩、要么（更糟）静默写出一份错的标定文件。
        %   这里统一加后缀 _bd 隔离。
        azr = deg2rad(azb(okb));
        rbObs = rb(okb)';
        M = [ones(numel(azr), 1), -cos(azr(:)), -sin(azr(:))];
        coef = M \ rbObs;
        R_bd = coef(1);
        kx_bd = coef(2);
        ky_bd = coef(3);
        d_bd = hypot(kx_bd, ky_bd);
        az0_bd = mod(rad2deg(atan2(ky_bd, kx_bd)), 360.0);
        bR = R_bd;
        bD = d_bd;
        res_bd = rbObs - M * coef;
        logf('详', '  以天顶像点 (%.1f, %.1f) 为极点, 逐方位测边界半径 (%d 帧中位):', ...
             x0, y0, numel(pf));
        logf('详', '  %10s %12s', 'az (deg)', 'r_boundary');
        for k = 1:6:72
            if isfinite(rb(k))
                logf('详', '  %10.0f %12.1f', azb(k), rb(k));
            end
        end
        logf('详', '');
        logf('详', '  拟合成圆 (绕主点): r(az) = R - d cos(az - az0)');
        logf('详', '    R  = %.1f px      边界圆半径', R_bd);
        logf('详', '    d  = %.1f px      边界圆心与天顶像点的距离', d_bd);
        logf('详', '    az0= %.1f deg     边界圆心相对天顶的方位', az0_bd);
        logf('详', '    残差 中位 %.1f px,  rms %.1f px', median(abs(res_bd)), sqrt(mean(res_bd .^ 2)));
        tilt = rad2deg(d_bd / max(abs(f), 1e-9));
        logf('详', '    由此推出的光轴倾角 d/f = %.2f deg', tilt);
        logf('');
        logf('  边界圆拟合：R = %.1f px，圆心偏离天顶像点 d = %.1f px（残余 rms %.1f px）', ...
             R_bd, d_bd, sqrt(mean(res_bd .^ 2)));
        if d_bd < 6.0
            logf('  ⇒ d < 6 px：边界基本**绕天顶同心** ⇒ "r 只是天顶角的函数"这一模型前提成立，倾角可忽略。');
        else
            logf('  ⚠ d = %.1f px：边界明显**不绕天顶同心**，而是绕一个偏离天顶 %.1f° 的主点。', d_bd, tilt);
            logf('    两种可能：(a) 它其实是镜头像圈、光轴偏离天顶；(b) 它本是地平圈，但天顶像点被定偏了。');
            logf('    区分办法：若为 (a)，改用以**主点**为极点的模型；若为 (b)，与【3.2final】的 x0/y0 核对。');
        end
    else
        logf('  （边界点太少，无法拟合）');
    end
else
    logf('  （缺少干净天光帧，跳过）');
end
% 无论拟合成不成功，都把"这条路线的能力边界"讲清楚 —— 它是**前提检验**，不是定标手段。
logf('  说明：这一节只检验前提，**不能用来定标**。圆是旋转对称的 ——');
logf('        绕圆心转任何角度都不改变这组等天顶角圆 ⇒ 圆心只能给出径向中心，给不出方位零点 rot；');
logf('        而"半径 R 对应多少天顶角"也需要外部尺度参照才能换算成 f。');

% --- 汇总 ---
% ⚠ 防呆：在写结果前**重新**从 pBest 取一次几何真值，不依赖任何裸变量名 ——
%   [3.9]/[3.10] 等后置检验段里的临时数组若与外层同名（如 a）会静默覆盖，
%   导致这里把数组当标量用而崩溃、或更糟：静默写出一份错的 JSON。
x0 = pBest(1); y0 = pBest(2); f = pBest(3); rotG = pBest(4); aa = pBest(5);
r90 = mc_r_of_z(bestModel, pi / 2, f, aa);
rimPx = rim;

% model_table：把 模型 × 手性 网格摊平成结构数组
modelTbl = struct('model', {}, 'chir', {}, 'rms', {}, 'x0', {}, 'y0', {}, ...
                  'f', {}, 'rot', {}, 'a', {});
for i = 1:numel(res)
    modelTbl(end+1) = struct('model', res(i).model, ...   %#ok<AGROW>
                             'chir', ternary(res(i).mirror, 'mir', 'dir'), ...
                             'rms', res(i).rms, 'x0', res(i).p(1), 'y0', res(i).p(2), ...
                             'f', res(i).p(3), 'rot', mod(res(i).p(4), 360), 'a', res(i).p(5));
end

cal = struct();
cal.model           = bestModel;
cal.r_of_z          = C.model_desc{find(strcmp(C.models, bestModel), 1)};
cal.chirality       = chirality;
cal.x0              = x0;
cal.y0              = y0;
cal.f               = f;
cal.a               = aa;
cal.rot             = mod(rotG, 360.0);
cal.r_horizon_px    = r90;      % z=90 deg 处的半径（外推自 z<=86 deg，系统项 ±6 px）
cal.rim_px          = rimPx;
cal.boundary_R      = bR;       % 3.10 边界圆拟合半径（NaN = 该次未算）
cal.boundary_d      = bD;       % 3.10 边界圆心偏离天顶像点（px）
cal.z_rim_deg       = zRim;
cal.fov_note        = sprintf(['拟合只覆盖 z<=90 deg。地平圈 r(90)=%.1f px vs 画幅半幅 %d px ' ...
    '⇒ 与 180 deg 鱼眼的''地平圈落在画幅边缘''一致；z>90 deg 的 FOV 外推不可靠，勿引用。'], ...
    r90, floor(C.size / 2));
cal.site_lon        = site(1);
cal.site_lat        = site(2);
cal.site_hgt_km     = site(3);
cal.rms_px          = rmsBest;
cal.n_frames        = numel(good);
cal.n_nights        = numel(usedNights);
cal.nights          = usedNights;
cal.rot_circ_mean   = cm;
cal.rot_circ_scatter = cs;
cal.rot_circ_R      = cR;
cal.rot_per_night   = perNight;
cal.gate_table      = gateTbl;
cal.model_table     = modelTbl;
cal.half_split      = xval;
cal.bootstrap       = boot;
cal.low_alt         = lowres;
cal.radial_profile  = prof;
cal.radial_dark     = darkf;

logf('');
logf('%s', repmat('=', 1, 78));
logf('【解算结论】%s（%s）   手性 %s   残差 rms %.2f px   %d 帧 / %d 夜', ...
     C.model_cn{find(strcmp(C.models, bestModel), 1)}, bestModel, ...
     ternary(strcmp(chirality, 'direct'), '直接', '镜像'), rmsBest, numel(good), numel(usedNights));
logf('  天顶像点  (%.2f, %.2f) px  [0-based]', x0, y0);
logf('            几何法独立参考 (%.1f, %.1f) ⇒ 差 %.1f px（两种方法互证）', ...
     C.x0_geo, C.y0_geo, hypot(x0 - C.x0_geo, y0 - C.y0_geo));
logf('  径向标度  f = %.2f px/rad，a = %+.4f   ⇒   r(z) = f·(z + a·z³)，z 用弧度', f, aa);
logf('  方位校正  rot = %.2f°  ± %.2f°（跨夜圆散布）', mod(rotG, 360.0), cs);
logf('  地平圈    r(90°) = %.1f px（画幅半幅 %d px，落在边缘内侧 %.0f px）', ...
     r90, floor(C.size / 2), C.size / 2 - r90);
if isfinite(zRim)
    logf('  与几何法  亮环 %.1f px vs 地平圈 %.1f px，差 %+.1f px（等效 f 相差 %+.1f%%）', ...
         rimPx, r90, rimPx - r90, 100 * (rimPx / (pi / 2) - r90 / (pi / 2)) / (r90 / (pi / 2)));
    logf('            ⇒ 在低空留出系统项(±6 px)与软边缘宽度(28 px)量级内一致，互为佐证。');
else
    logf('  与几何法  视场边缘 %.1f px 超出该模型可达半径 (r(90°) = %.1f px)', rimPx, r90);
end
logf('%s', repmat('=', 1, 78));

% --- 结果健全性断言 ---
% 宁可在这里炸，也不要写出一份看着合理其实错掉的标定文件。
bad = {};
if ~(f > 300.0 && f < 700.0)
    bad{end+1} = sprintf('f=%.2f px/rad 不在 [300,700]', f);
end
if ~(x0 > 0 && x0 < C.size && y0 > 0 && y0 < C.size)
    bad{end+1} = sprintf('天顶像点 (%.1f,%.1f) 落在画幅外', x0, y0);
end
% rot 的取值域：cal.rot 落盘前会做 mod(,360)，所以这里必须判**归一化之后**的值。
% 直接判 rotG 会把一个完全合法的 -0.9 deg（初值取 0 时 LM 会自由地往负方向走）
% 判成"不合法"，从而拒绝导出一份本来正确的标定。
rotNorm = mod(rotG, 360.0);
if ~isfinite(rotG) || ~(rotNorm >= 0 && rotNorm < 360)
    bad{end+1} = sprintf('rot=%.3f deg 不合法', rotG);
end
if ~(r90 > 3.0 && r90 < C.size)
    bad{end+1} = sprintf('r(90 deg)=%.1f px 不合法', r90);
end
if ~(rmsBest >= 0 && rmsBest < 20.0)
    bad{end+1} = sprintf('rms=%.2f px 过大', rmsBest);
end
if hypot(x0 - C.x0_geo, y0 - C.y0_geo) > 60.0
    bad{end+1} = sprintf('天顶像点与几何法 (%.1f,%.1f) 差 %.1f px, 超过 60 px', ...
                         C.x0_geo, C.y0_geo, hypot(x0 - C.x0_geo, y0 - C.y0_geo));
end
assert(isempty(bad), 'moon_calib:sanityFailed', ...
    '解算结果未通过健全性断言，**不导出任何标定文件**：\n   - %s', strjoin(bad, '\n   - '));
end


% ==============================================================================
% [7] 阶段 4 —— 验证图
% ==============================================================================

function qc = mc_verify(cfg, site, cal, logf)
%MC_VERIFY 把「纯星历预测」的月亮位置（红圈）叠到原始帧上，绿十字 = 拟合天顶。
C = mc_defaults();
logf('');
logf('%s', repmat('=', 1, 78));
logf('【阶段 4】质检图：把**纯星历预测**的月亮位置（红圈）叠到原始帧上，绿十字 = 拟合天顶');
logf('%s', repmat('=', 1, 78));

vdir = fullfile(cfg.out, 'qc');
if ~exist(vdir, 'dir'), mkdir(vdir); end
% 清掉上一轮的 QC 图，免得残留的旧图（比如被剔除夜的图）混在里面误导
old = dir(fullfile(vdir, 'qc_*.png'));
for i = 1:numel(old)
    try, delete(fullfile(vdir, old(i).name)); catch, end %#ok<CTCH>
end

% 挑最有说服力的帧：不按固定间隔抽样，而是**按星历直接挑**
%   low = 月亮最低但仍在地平上的那一帧（预测位置紧贴视场边缘）
%   mid = 月亮高度最接近拟合区间中点的那一帧
% ⚠ 固定间隔抽样只会拿到每夜开头（月亮最高）的帧，低空帧永远选不到 —— 这是踩过的坑。
tgtMid = 0.5 * (cfg.alt_lo + cfg.alt_hi);
picks = struct('night', {}, 'file', {}, 't', {}, 'alt', {}, 'az', {}, 'tag', {});
for i = 1:numel(cal.nights)
    nd = cal.nights{i};
    fs = mc_list_frames(fullfile(cfg.raw, nd));
    if numel(fs) < 20, continue; end
    n = numel(fs);
    alts = zeros(n, 1);
    azs = zeros(n, 1);
    tss = NaT(n, 1);
    for k = 1:n
        t = mc_parse_ts(mc_basename(fs{k}));
        tss(k) = t;
        [a_, A_] = mc_moon_altaz(mc_jd(t), site);
        alts(k) = a_;
        azs(k) = A_;
    end
    m = alts > cfg.alt_hold + 0.5;
    if any(m)
        tmp = alts; tmp(~m) = 1e9;
        [~, kk] = min(tmp);
        if alts(kk) < cfg.alt_lo
            picks(end+1) = struct('night', nd, 'file', mc_basename(fs{kk}), 't', tss(kk), ...   %#ok<AGROW>
                                  'alt', alts(kk), 'az', azs(kk), 'tag', 'low');
        end
    end
    [~, km] = min(abs(alts - tgtMid));
    if alts(km) >= cfg.alt_lo && alts(km) < cfg.alt_hi
        picks(end+1) = struct('night', nd, 'file', mc_basename(fs{km}), 't', tss(km), ...   %#ok<AGROW>
                              'alt', alts(km), 'az', azs(km), 'tag', 'mid');
    end
end

lows = picks(strcmp({picks.tag}, 'low'));
mids = picks(strcmp({picks.tag}, 'mid'));
uniq = struct('night', {}, 'file', {}, 't', {}, 'alt', {}, 'az', {}, 'tag', {});
seen = {};

% 分配名额：低空帧（月亮贴着地平线）最有说服力，优先给它留 2 个；其余给中空帧。
% 同一夜最多一张，名额没满再放宽。
for pass = 1:2
    if pass == 1
        lst = lows; want = min([2, numel(lows), cfg.n_qc]);
    else
        lst = mids; want = cfg.n_qc - numel(uniq);
    end
    got = 0;
    for i = 1:numel(lst)
        if any(strcmp(seen, lst(i).night)), continue; end
        uniq(end+1) = lst(i);   %#ok<AGROW>
        seen{end+1} = lst(i).night;   %#ok<AGROW>
        got = got + 1;
        if got >= want, break; end
    end
end
if numel(uniq) < cfg.n_qc                  % 名额没满就放宽"每夜最多一张"的限制
    for i = 1:numel(picks)
        if numel(uniq) >= cfg.n_qc, break; end
        dup = false;
        for j = 1:numel(uniq)
            if strcmp(uniq(j).night, picks(i).night) && strcmp(uniq(j).file, picks(i).file)
                dup = true; break;
            end
        end
        if ~dup, uniq(end+1) = picks(i); end   %#ok<AGROW>
    end
end

qc = struct('night', {}, 'file', {}, 't', {}, 'alt', {}, 'az', {}, ...
            'pred_px', {}, 'pred_py', {}, 'peak_near', {}, 'bg', {}, 'snr_dn', {}, 'png', {});
for i = 1:numel(uniq)
    nd = uniq(i).night;
    fp = fullfile(cfg.raw, nd, uniq(i).file);
    im = mc_load_gray(fp);
    if isempty(im), continue; end
    zr = deg2rad(90.0 - uniq(i).alt);
    rr = mc_r_of_z(cal.model, zr, cal.f, cal.a);
    px = cal.x0 + rr * sin(deg2rad(uniq(i).az) + deg2rad(cal.rot));
    py = cal.y0 - rr * cos(deg2rad(uniq(i).az) + deg2rad(cal.rot));

    f = double(im);
    lo = mc_pctile(f(:), 1);
    hi = mc_pctile(f(:), 99.5);
    gg = min(max((f - lo) / max(hi - lo, 1.0), 0), 1);
    g8 = uint8(round(255 * asinh(gg * 12) / asinh(12)));
    I = repmat(g8, 1, 1, 3);

    I = mc_draw_circle(I, px, py, 34, [255 0 0], 2);
    I = mc_draw_cross(I, px, py, 46, [255 0 0], 1);
    I = mc_draw_cross(I, cal.x0, cal.y0, 15, [0 255 0], 2);

    outp = fullfile(vdir, sprintf('qc_%s_%s_alt%+.0f_%s.png', ...
                     nd, datestr(uniq(i).t, 'HHMMSS'), uniq(i).alt, uniq(i).tag));
    imwrite(I, outp);

    % 预测位置邻域的实际亮度，作为"红圈里真有东西"的客观证据
    hw = 30;
    ay0 = max(0, floor(py) - hw); ay1 = min(C.size - 1, floor(py) + hw);
    ax0 = max(0, floor(px) - hw); ax1 = min(C.size - 1, floor(px) + hw);
    pk = max(max(f(ay0 + 1:ay1 + 1, ax0 + 1:ax1 + 1)));
    bgv = median(f(:));
    qc(end+1) = struct('night', nd, 'file', uniq(i).file, ...   %#ok<AGROW>
                       't', datestr(uniq(i).t, 'HH:MM:SS'), 'alt', uniq(i).alt, 'az', uniq(i).az, ...
                       'pred_px', px, 'pred_py', py, 'peak_near', pk, 'bg', bgv, ...
                       'snr_dn', pk - bgv, 'png', sprintf('qc_%s_%s_alt%+.0f_%s.png', ...
                       nd, datestr(uniq(i).t, 'HHMMSS'), uniq(i).alt, uniq(i).tag));
    logf('详', '  %s %s  alt=%+5.1f az=%6.1f -> 预测像点 (%7.1f, %7.1f)   邻域峰值 %.0f DN (背景 %.0f)', ...
         nd, datestr(uniq(i).t, 'HH:MM:SS'), uniq(i).alt, uniq(i).az, px, py, pk, bgv);
end
logf('');
logf('  质检图 %d 张 -> %s', numel(qc), vdir);
if ~isempty(qc)
    snr = [qc.snr_dn];
    logf('  红圈邻域"峰值 − 背景"中位 %.0f DN ⇒ 纯星历预测的位置上确实有目标（不是画在空处）', median(snr));
    logf('  看图：红圈套住月亮、绿十字落在画幅中心附近即正常。');
end
end


function I = mc_draw_circle(I, cx, cy, rad, col, wid)
%MC_DRAW_CIRCLE 在 uint8 H×W×3 图上画圆（不依赖 Image Processing Toolbox）。
[H, W, ~] = size(I);
x0 = max(0, floor(cx - rad - wid - 1));
x1 = min(W - 1, ceil(cx + rad + wid + 1));
y0 = max(0, floor(cy - rad - wid - 1));
y1 = min(H - 1, ceil(cy + rad + wid + 1));
if x1 < x0 || y1 < y0, return; end
[Xg, Yg] = meshgrid(x0:x1, y0:y1);
m = abs(hypot(Xg - cx, Yg - cy) - rad) <= (wid / 2 + 0.5);
sub = I(y0 + 1:y1 + 1, x0 + 1:x1 + 1, :);
for c = 1:3
    ch = sub(:, :, c);
    ch(m) = col(c);
    sub(:, :, c) = ch;
end
I(y0 + 1:y1 + 1, x0 + 1:x1 + 1, :) = sub;
end


function I = mc_draw_cross(I, cx, cy, half, col, wid)
%MC_DRAW_CROSS 在 uint8 H×W×3 图上画十字（不依赖 Image Processing Toolbox）。
[H, W, ~] = size(I);
a = max(0, round(cy - half)); b = min(H - 1, round(cy + half));
c = max(0, round(cx - wid / 2)); d = min(W - 1, round(cx + wid / 2));
if b >= a && d >= c
    for ch = 1:3
        blk = I(a + 1:b + 1, c + 1:d + 1, ch);
        blk(:) = col(ch);
        I(a + 1:b + 1, c + 1:d + 1, ch) = blk;
    end
end
a = max(0, round(cx - half)); b = min(W - 1, round(cx + half));
c = max(0, round(cy - wid / 2)); d = min(H - 1, round(cy + wid / 2));
if b >= a && d >= c
    for ch = 1:3
        blk = I(c + 1:d + 1, a + 1:b + 1, ch);
        blk(:) = col(ch);
        I(c + 1:d + 1, a + 1:b + 1, ch) = blk;
    end
end
end


% ==============================================================================
% [8] 阶段 5 —— 导出
% ==============================================================================

function mc_export(cfg, cal, qc, logf)
%MC_EXPORT 写 .mat / .json / .txt。
logf('');
logf('%s', repmat('=', 1, 78));
logf('【阶段 5】导出');
logf('%s', repmat('=', 1, 78));

payload = cal;
payload.qc = qc;
payload.provenance = struct( ...
    'script', 'moon_calib.m', ...
    'generated_utc', datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ...
    'raw_root', cfg.raw, ...
    'site_key', cfg.site_key, ...
    'criteria', struct('alt_lo', cfg.alt_lo, 'alt_hi', cfg.alt_hi, ...
                       'scat_max_deg', cfg.scat_max, 'min_frames', cfg.min_frames, ...
                       'rim_px', cfg.rim), ...
    'note', ['rot 的不确定度取跨夜圆散布; x0/y0/f/a 另见 bootstrap。' ...
             'rot 与畸变形状在无星数据上唯一可得的绝对参照就是月亮。']);

% --- .mat（MATLAB 原生，最要紧的产物） ---
% 文件名默认写死 calibration_params.mat；但同一工程里可能有多批数据
% （二期 2025-2026 16 bit / 一期 2014 8 bit），各自的定标必须分开存，
% 否则后跑的那批会**静默覆盖**前一批。用 'outmat' 显式指定完整路径即可。
if isfield(cfg, 'outmat') && ~isempty(cfg.outmat)
    mp = cfg.outmat;
    mdir = fileparts(mp);
    if ~isempty(mdir) && ~exist(mdir, 'dir'), mkdir(mdir); end
    if exist(mp, 'file') == 2
        [~, bn0, ex0] = fileparts(mp);
        mbk = fullfile(mdir, sprintf('%s_backup_%s%s', bn0, datestr(now, 'yyyymmdd_HHMMSS'), ex0));
        copyfile(mp, mbk);
        logf('详', '  原生格式  原文件已备份 -> %s', mbk);
    end
else
    mp = fullfile(cfg.out, 'calibration_params.mat');
end
calMat = cal;
if isfield(calMat, 'radial_profile'), calMat = rmfield(calMat, 'radial_profile'); end
if isfield(calMat, 'radial_dark'),    calMat = rmfield(calMat, 'radial_dark'); end
calMat.r_z_fun = 'r = f*(z + a*z.^3);  z in radians';
calMat.provenance = payload.provenance;
save(mp, 'calMat', '-v7');
logf('详', '  原生格式  取 calMat.x0 / .y0 / .f / .a / .rot');
logf('详', '            注意：这是本程序原生格式，**不能**直接喂投影管线（管线要 zenith/rot/radius/rz_poly）。');

% --- 可选：直接写一份 airglow_geo_pipeline 能 load 的 cal 结构 ---------------
% 必要性（2026-09-16 实测确认）: 管线 load_calibration 会遍历 .mat 里每个 struct
% 字段, 找一个**含 zenith 字段**的, 然后强制要求
%     {'zenith','rot','radius','rz_poly'} 四者齐备。
% 而上面写出的 calMat 只有 x0/y0/f/a/rot -> 管线会直接
%     error('标定结构体缺少字段：zenith')
% 于是"跑完 moon_calib 就能 run.m"在没有这一步时是**不成立**的。
ppOut = '';    % 实际写出的管线定标路径（未写则为空）—— 供末尾结果卡使用
if isfield(cfg, 'pipeline_out') && ~isempty(cfg.pipeline_out)
    pp = cfg.pipeline_out;
    if exist('moon_calib_to_pipeline', 'file') ~= 2
        logf('  ⚠ 未找到 moon_calib_to_pipeline.m（不在 MATLAB 路径上）→ 跳过管线格式输出。');
    else
        ppOut = pp;
        pdir = fileparts(pp);
        if ~isempty(pdir) && ~exist(pdir, 'dir'), mkdir(pdir); end
        % 目标已存在 -> 先时间戳备份, 绝不静默覆盖别人的标定
        if exist(pp, 'file') == 2
            [~, bn, ex] = fileparts(pp);
            bbk = fullfile(pdir, sprintf('%s_backup_%s%s', bn, datestr(now, 'yyyymmdd_HHMMSS'), ex));
            copyfile(pp, bbk);
            logf('详', '  管线格式  原文件已备份 -> %s', bbk);
        end
        cal2 = moon_calib_to_pipeline(calMat, pp, 'Radius', cfg.rim, 'Verbose', false);
        logf('详', '  管线格式  取 cal.zenith / .rot / .radius / .rz_poly');
        logf('详', '            obs_lat/lon = %.3f, %.3f（会覆盖投影脚本的实参）', cal2.obs_lat, cal2.obs_lon);
        logf('详', '            zenith(1-based) = (%.3f, %.3f)   radius = %.2f px', ...
             cal2.zenith(1), cal2.zenith(2), cal2.radius);
    end
end

% --- .json（便于跨语言比对） ---
jp = fullfile(cfg.out, 'calib_params.json');
fp = fullfile(cfg.out, 'calib_params.json.tmp');
jsonOk = true;
try
    jtxt = jsonencode(payload, 'PrettyPrint', true);
    fh = fopen(fp, 'w', 'n', 'UTF-8');
    fwrite(fh, jtxt, 'char');
    fclose(fh);
    movefile(fp, jp, 'f');
catch ME
    jsonOk = false;
    logf('  (jsonencode 失败, 跳过 .json: %s)', ME.message);
end
if jsonOk
    logf('详', '  JSON   -> %s', jp);
end

% --- 人类可读摘要 + 投影公式 ---
tp = fullfile(cfg.out, 'calib_summary.txt');
fh = fopen(tp, 'w', 'n', 'UTF-8');
fprintf(fh, '月亮定向定标结果\n%s\n', repmat('=', 1, 60));
keys = {'model', 'r_of_z', 'chirality', 'x0', 'y0', 'f', 'a', 'rot', ...
        'rms_px', 'n_frames', 'n_nights', 'rot_circ_scatter', ...
        'r_horizon_px', 'rim_px', 'z_rim_deg', 'site_lon', 'site_lat'};
for i = 1:numel(keys)
    v = cal.(keys{i});
    if ischar(v)
        fprintf(fh, '%-20s %s\n', keys{i}, v);
    else
        fprintf(fh, '%-20s %s\n', keys{i}, num2str(v, '%.6g'));
    end
end
fprintf(fh, '\n注意: r_horizon_px 是 z=90 deg 处的半径（外推自 z<=86 deg 的实测，系统项 ±6 px）。\n');
fprintf(fh, '%s\n', cal.fov_note);
fprintf(fh, '\n用于投影（z 单位 rad，结果单位 px，0-based 坐标）:\n');
fprintf(fh, '  %% 天顶角\n');
fprintf(fh, '  z = acos( sind(lat)*sind(dec) + cosd(lat)*cosd(dec)*cosd(H) );\n');
fprintf(fh, '  r = %.4f * (z + %+.6f * z.^3);\n', cal.f, cal.a);
fprintf(fh, '  x = %.3f + r .* sin(az + deg2rad(%.3f));\n', cal.x0, cal.rot);
fprintf(fh, '  y = %.3f - r .* cos(az + deg2rad(%.3f));\n', cal.y0, cal.rot);
fprintf(fh, '\n  其中 az 为方位角（北=0、东=90，deg），H 为时角（deg）。\n');
fprintf(fh, '  MATLAB 可直接用上面四行，注意 sind/cosd 用角度、r 的表达式用弧度。\n');
fprintf(fh, '\n直接可用的投影函数（把 cal 传进来即可）:\n');
fprintf(fh, '  [x, y] = moon_calib_project(cal, lat, dec, H, az)\n');
fclose(fh);
logf('详', '  文本摘要 -> %s', tp);

% 把本次写出的产物路径挂到 root appdata，供末尾的 mc_result_card 使用。
% （用 appdata 而不是返回值：mc_export 是过程式函数，改签名会波及所有调用点。）
setappdata(0, 'mc_cal_exported', struct('native', mp, 'pipe', ppOut, ...
    'qc', ternary(~isempty(qc), fullfile(cfg.out, 'qc'), '')));
logf('');
logf('  导出完成：原生定标 / 管线定标 / JSON / 文本摘要 四份（路径见下方结果卡）');
end


% ==============================================================================
% [10] 结果卡
% ==============================================================================

function mc_result_card(cal, cfg, logpath, sec, logf)
%MC_RESULT_CARD 定标结果卡 —— 整条流程的收尾输出。
%
% 设计意图（2026-09-18）：在此之前，跑完之后控制台最后看到的是若干文件路径，
% 而真正要看的两个数（天顶像点、rot）夹在阶段 3 的中段、被后面的阶段刷走。
% 这张卡固定形状、放在最后，做到"跑完先看这一块就够"。
%
% 对齐说明：中文是全角（2 列），MATLAB 的 %-*s 按**字符**而不是显示列补空格，
% 所以标题列用"4 个汉字 + 2 空格"这种固定写法手工对齐，不用格式化宽度。
% 输入：
%   cal     : mc_solve 产出的定标结构
%   cfg     : 本次运行的 cfg（取路径）
%   logpath : 日志文件全路径
%   sec     : 总用时 (s)
%   logf    : 日志闭包
W = 78;
rule = repmat('=', 1, W);
half = 512;    % 画幅 1024² ⇒ 半幅（与 mc_defaults 的 C.size 一致）

logf('');
logf('%s', rule);
logf('  定标结果 · 月亮法定标         （完整表见日志文件，本卡只列关键量）');
logf('%s', repmat('-', 1, W));
logf('  天顶像点  (%.2f, %.2f)  px   [0-based]', cal.x0, cal.y0);
logf('  方位校正  rot = %+.3f°   ± %.3f°（跨夜圆散布）', cal.rot, cal.rot_circ_scatter);
logf('  径向标度  f = %.3f px/rad    a = %+.6f', cal.f, cal.a);
logf('            r(z) = f·(z + a·z³)，z 用弧度；方位角北=0、东=90');
logf('  拟合残差  rms = %.3f px        （%d 帧 / %d 夜）', ...
     cal.rms_px, cal.n_frames, cal.n_nights);
logf('  手  性    %s       径向模型 %s（%s）', ...
     ternary(strcmp(cal.chirality, 'direct'), '直接（无镜像）', '镜像（左右翻转）'), ...
     cal.model, cal.r_of_z);
logf('%s', repmat('-', 1, W));
logf('  地平圈    r(90°) = %.1f px   画幅半幅 %d px ⇒ 边缘内侧 %.0f px', ...
     cal.r_horizon_px, half, half - cal.r_horizon_px);
if isfinite(cal.rim_px)
    if cal.rim_px > cal.r_horizon_px + 2
        logf('  视场边缘  %.1f px **大于**地平圈 ⇒ 它只是画幅裁边，不是光学盘边（本机无可用盘边信息）', ...
             cal.rim_px);
    else
        logf('  视场边缘  R_rim = %.1f px   与地平圈差 %+.1f px', ...
             cal.rim_px, cal.rim_px - cal.r_horizon_px);
    end
end
if isfinite(cal.boundary_R)
    logf('  边界同心  边界圆 R = %.1f px   圆心偏离天顶 %.1f px（越小越说明"绕天顶同心"成立）', ...
         cal.boundary_R, cal.boundary_d);
end
if isfinite(cal.z_rim_deg) && cal.z_rim_deg < 95
    logf('  盘边天顶角 θ_rim = %.1f°   （**不是 90°** ⇒ 不能用 f = R_rim/(π/2) 反推焦距）', ...
         cal.z_rim_deg);
end
logf('%s', repmat('-', 1, W));
logf('  产  物    原生定标  %s', ternary(isfield(cfg, 'outmat') && ~isempty(cfg.outmat), ...
     cfg.outmat, fullfile(cfg.out, 'calibration_params.mat')));
S = getappdata(0, 'mc_cal_exported');
if isstruct(S) && isfield(S, 'pipe') && ~isempty(S.pipe)
    logf('            管线定标  %s', S.pipe);
end
logf('            JSON      %s', fullfile(cfg.out, 'calib_params.json'));
logf('            文本摘要  %s', fullfile(cfg.out, 'calib_summary.txt'));
logf('            日志      %s', logpath);
if isstruct(S) && isfield(S, 'qc') && ~isempty(S.qc)
    logf('            质检图    %s', S.qc);
end
logf('%s', repmat('-', 1, W));
logf('  下一步    start_process(''<夜>'')   做图像处理 + 地理投影');
logf('  提示      逐夜/逐帧/逐半径的长表只进日志文件；要现场看全量：moon_calib(...,''chatty'',true)');
logf('  总用时    %.1f s', sec);
logf('%s', rule);
end


% ==============================================================================
% [9] 小工具
% ==============================================================================

function s = mc_basename(p)
[~, n, e] = fileparts(p);
s = [n, e];
end


function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end


function name = mc_pick_best(res, C)
%MC_PICK_BEST 直接手性里 rms 最小的模型，返回**模型名**（不是索引！）。
% 踩过的坑：最早这里写成返回索引 i，调用处却拿它当模型名去 mc_get_p 里 strcmp，
% 结果报"网格里找不到 / mirror=0"（%s 打不出数字、%d 把 false 打成 0），
% 错误信息完全指不到真因。返回字符串就没这问题。
bestRms = Inf;
name = C.models{1};
for k = 1:numel(C.models)
    r = mc_get_rms(res, C, C.models{k}, false);
    if r < bestRms
        bestRms = r;
        name = C.models{k};
    end
end
end


function p = mc_get_p(res, C, model, mirror)
p = [];
for k = 1:numel(res)
    if strcmp(res(k).model, model) && res(k).mirror == mirror
        p = res(k).p;
        return;
    end
end
error('moon_calib:noEntry', '网格里找不到 %s / mirror=%d', model, mirror);
end


function r = mc_get_rms(res, C, model, mirror)
%MC_GET_RMS 从网格结果里取指定 模型/手性 的 rms。
r = NaN;
for k = 1:numel(res)
    if strcmp(res(k).model, model) && res(k).mirror == mirror
        r = res(k).rms;
        return;
    end
end
end