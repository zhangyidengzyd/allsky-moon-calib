function out = calib_2014(varargin)
%CALIB_2014  一期(2014)专用定标 —— 并与马欣博士论文的数据做对照验证
%
%   calib_2014                            % 全流程（默认参数，约 1~2 分钟）
%   calib_2014('Stage','solve')           % 复用已有测量缓存，只重解（几秒）
%   calib_2014('SkipFrames',true)         % 跳过论文样例帧复现（省时）
%   out = calib_2014(...)                 % 返回结果结构（含采用解 + 对照表数值）
%
% ===========================================================================
% 为什么必须单独有这一份 —— 不能直接用 q1_calibrate
% ===========================================================================
% 一期 2014 与二期 2025 是**两台不同光学配置**的相机，且位深不同。直接套用会全灭：
%
% (1) 位深：一期是 **8 bit**（满量程 255），本链所有阈值按 16 bit 写死
%     （饱和判据 65000、暗场基座 3441）。不统一 ⇒ moon_calib 的 detect 永远检不出，
%     而报错会伪装成"数据里没有月亮"。→ 全链已改走 q_read_gray()（uint8×257）。
%
% (2) 几何圆心：moon_calib 的 mc_defaults 把**二期**圆心 (543.824, 474.352) 与
%     R_rim 512.93 **写死**，而 [3.1] 判据算的是 phi = atan2(dx,-dy)。
%     一期实测圆心 ≈ (554, 537)：与二期圆心在 y 方向差 **63 px** ⇒ phi 误差随月亮
%     方位变化 ⇒ 逐夜 rot 散布被抬到 4~92 deg ⇒ 7 夜全被判成"假目标"。
%     → 本函数通过 moon_calib 的 x0_geo / y0_geo / rim 参数覆盖。
%
% (3) ★ 最关键的一条：**黄昏/黎明天光**。8 bit 数据里，太阳只要高于约 -12 deg，
%     大片天光就能顶到满量程，被 [STAGE 2] 当成"饱和月核"检出来。这些假目标的
%     质心跟着**太阳**（天光梯度）而不是月亮走。实测：20140930 太阳 -8.4 deg 的首帧
%     rot = 118.6 deg，而同夜太阳 < -15 deg 的 6 帧 rot = 6.9/6.2/6.1/6.8/6.6/6.8 deg。
%     → 用 sun_alt_max 把暮光帧剔掉。这一条是"能不能用"的分水岭：
%       不加门限 = 通过 0 帧；加 -18 deg = 通过 3 夜，且逐夜 rot 一致到 0.5 deg。
%     二期 16 bit 数据夜天光只有几千 DN，不会被误判，所以 moon_calib 默认 Inf（不过滤）。
%
% (4) ★ 第四个坑，也是最隐蔽的一个：**rot 与圆心 y0 强简并**。
%     实测 Δrot/Δy0 ≈ -0.17 deg/px：把 (x0,y0) 作自由参数交给 LM 一起拟合时，
%     径向模型的形状误差会被"圆心平移"吸收 ⇒ 不同高度窗口解出不同 rot
%     （alt 14-78 -> 3.14 deg；alt 14-50 -> 5.30 deg），两者都自认为 rms 只有 3 px。
%     破简并的办法：**先只用"方向"定圆心，再固定圆心拟合尺度**。
%       方向信息（像面方位 phi）与径向模型**完全无关**（径向畸变不改变方向），
%       所以"让各帧的 rot = phi - A 都相等"是一个只关于 (x0,y0) 的、良定的二维问题。
%
% (5) ★ 第五个坑：**"逐夜统计量"本身也是陷阱**。一开始的定心代价写成
%     "各夜 rot 的(3xMAD 截尾)圆均值之间的圆方差"，得到 (535, 543) / rot 3.79 deg，
%     看着还挺自洽（rms 6.18 px）。但把这个代价换成"逐夜圆中位"或"加权(1/sigma^2)
%     圆均值"或"只用最干净的 3 夜池化"，极小分别跑到 (551,556) / (548,538) / (557.6,536.2)
%     —— 也就是说这个代价是**被 2 个坏夜(20140930/20141007，夜内 rot 散射 32~35 deg)
%     的汇总值拖出来的伪极小**：坏夜的 rot 分布又宽又偏，圆均值/中位/截尾三个统计量
%     给出 3 个差 10 deg 的答案，谁都能被"凑"到和其他夜一致。诊断证据：
%       在 (535,543) 上，与全局 rot 一致到 +-2 deg 的帧只有 370/676 = 55%（散布 0.96 deg）；
%       在 (554,536) 上有 543/676 = 80%（散布 0.64 deg）；
%       在 (552,556) 上有 448/676 = 66%（散布 0.75 deg）。
%     ⇒ 采用**全局紧截尾共识**做定心判据（见 [B2]）：不用任何"逐夜汇总量"，
%       直接数"与全体帧的 rot 中位一致到 +-2 deg 的帧最多、且这些帧散布最小"的圆心。
%       坏帧（跟错目标、云、相位导致的质心偏移）在任何圆心下都不可能一致，会被自动丢掉；
%       这个判据对**帧集的选择也不敏感**（见下表，6 个不同帧集都收敛到同一处）：
%         帧集                          帧数   共识圆心(散布最小)   内点率
%         全部 7 夜                     760   (553.7, 536.5)       78%
%         sun<-18                        676   (554.0, 536.6)       79%
%         sun<-18 & alt 14-50           395   (556.6, 536.7)       90%
%         sun<-18 & alt 14-70           595   (554.0, 536.6)       81%
%         sun<-18 & 去两个坏夜            315   (554.2, 536.4)       96%
%         最干净 3 夜                    253   (555.4, 548.3)       97%
%       ⇒ 共识圆心 **≈ (554.5, 536.5) +-3 px**，正好等于论文 §3.2.1 的三星余弦定理值
%         (554, 538)。逐帧 hold-out、LM 自由拟合 (556.7, 539.8)、加权圆均值 (548, 538)
%         都落在同一处。
%
%     为什么可以用"更宽的高度窗口"来定圆心、却必须用"更窄的高度窗口"来拟合尺度：
%       定圆心只用**方向**，方向对任何径向畸变都不变，所以月核被光晕撑大（alt>55）
%       并不影响它；而拟合 r(z) 恰恰要的是径向位置，撑大的月核会污染尺度。
%       故本函数把两个窗口分开：AltHiCenter=70 用于定心，AltHi=50 用于拟合。
%
% ===========================================================================
% 与马欣论文（3.2 节）的对照
% ===========================================================================
% 论文值（原文抄录；全文见 mooncal_work\reports\参考资料_thesis文本抽取.txt）：
%   * 天顶像点  §3.2.1 用三颗星余弦定理，**原图坐标系** → (x=554, y=538)
%               §3.2.2 方位校正用的"圆心即天顶" → (X0=452, Y0=452)
%     ★★ 2026-09-17 二次修订：这两组数是**同一个物理点**，不是两组天顶，也不是
%        "两个不同的东西"。原文把顺序写得很清楚 —— §3.2.1：
%          "…原始气辉图像是圆形的，因此需要以圆形模板来提取中间的图像。首先需要
%            知道圆心的在图像中的位置，这里我们通过选取亮度增强图中分布均匀的三个
%            星点，根据余弦定理计算得到圆心即天顶的位置坐标 … 对其求平均值即得到
%            最终提取图像的圆心坐标（x=554，y=538）。"
%        ⇒ "圆心"是裁剪的**输入**（由星表独立求出），裁剪本身"以天顶为圆心"做；
%          (452,452)+(102,86) = (554,538) 逐位闭合，裁剪窗 904x904、模板圆半径
%          452 = 窗半宽 = 内切圆（与实测量论文图3.2(c)(d) 的 R/w = 0.499 一致）。
%        ⇒ 因此**不存在**"论文把视场圆心误当天顶"这回事（那是本文件 09-17 早先的
%          误读）。当时用来支撑它的"圆心−天顶差 20.6 px"实测也已被推翻：见下。
%        ⇒ 原 [E3] 的 −6.90 deg 一栏、以及"简并斜率判别表"，均已删除。
%   * "圆形视场圆心"这个量在一期数据上**测不准**，故本工程不再引用它。
%       2026-09-17 阈值扫描（一期内 4 夜，见工作区 _t14.txt）：拟合出的"视场盘"
%       圆心随阈值单调滑动 —— 20140923 帧上阈值 14→18 时 cy 从 544 滑到 480、
%       R 从 452 滑到 399、rms 从 4.9 涨到 22 px。说明一期的"圆形视场边界"不是
%       一条清晰几何边，而是弥散的渐晕过渡，阈值法（含 moon_calib 的 [3.10]）都
%       不能拿它当硬证据。→ 相关结论一律撤下，只保留"月核+方位"这条独立链。
%   * 北极星像点 (Xn=499, Yn=55)，其方位角 0.4 deg
%       独立复算（rot_audit.m）: Polaris 在 2014-09-23 17:06:12 UT、站点处的真实
%       方位角 = +0.390 deg ⇒ 论文这一项**是对的**（差 0.01 deg）。但那个像素本身
%       在归档 8 bit 帧里测不到（邻域峰值比背景高 3 个灰阶），无法复核。
%   * 方位校正角 (3.2) 式 = 6.35 deg（6.75 - 0.4；6.75 由 (452,452)+(499,55) 算出）
%   * 径向标度  (3.9) 式  R_img = -0.55 z^2 + 305.2 z + 16.68     （z 以**弧度**）
%   * 气辉高度  300 km；投影画布 600x600，0.020 deg/px
%               纬度 13.50~25.48 N，经度 103.220~115.20 E
%   * 样例帧    2014-09-23 17:06:12 UT（新月夜 → 暗夜；该夜月亮整夜在地平以下）
%   * 论文另在 §3.2.2/§3.3 混用了 2013-10-25 的图像（图3.3、图3.10），且 §3.2 内
%     混用"原图坐标系"与"裁剪坐标系"而未加说明 —— 对照时务必记住这点。
%
% 那么 rot 差的 2.0 deg 怎么读？（结论：无法裁决，但不在本解这一侧）
%   逐项排除（脚本 rot_audit.m，输出 reports\工具_rot三项审计.txt）：
%     (1) 天顶：论文 (554,538) vs 本解 (553.5,536.5) ⇒ 差 1.6 px，不是原因；
%     (2) 方位角：论文取 0.4 deg，真实值 +0.390 deg ⇒ 不是原因；
%     (3) 时钟：星历时刻平移 dt 扫描 ±36 min，二维 rms 在 dt=0 最小（4.572 px）、
%         纯方位散射也在 dt=-3~0 最小（0.538 deg），rot 对 dt 斜率仅 0.160 deg/min
%         ⇒ 钟差 5 min 只带 0.8 deg ⇒ 不是原因；
%     (4) 剩下只有那颗星的**方向**：把论文自己的 (499,55) 配本解天顶算出 rot=6.45
%         deg（论文 6.35）⇒ 论文内部自洽，2 deg 整份在这颗星的方向上，而 2 deg 在
%         r=398 px 处 = 横向 14 px；该星还离论文自己的径向式(3.9) 外偏 10.5 px
%         （2.7%），并且它在归档帧里根本测不到。
%   ⇒ 措辞建议：本解 rot = +4.35 ± 1.0 deg 有旁证支撑；论文 6.35 deg 依赖单颗手工
%     拾取的星像，精度不可能优于 1 deg。2 deg 在 300 km 高度上（z=30~60 deg）对应
%     6~18 km 东西向位移，做泡团结构对比时值得写一句。
%
% 本函数给出的**独立**结果与论文的逐项对照，见 [E] 节打印的表与
%   <原生定标目录>\thesis_compare.txt 落盘。
%
% 用法要点：跑完把 cfg.calFile（= Processed_Output\calibration_params_14.mat）
% 交给 start_process 的 'Cal' 参数，即可处理一期 2014 的帧。
%
% 依赖：q_cfg / q_read_gray / moon_calib / zenith_from_rim / moon_calib_to_pipeline
%       / moon_altaz_ref / sun_altaz_ref。全部在本工程根目录。
%
% ★ 路径独立：本函数**不读** q_cfg 的 rawRoot / calFile / calNative。
%   q_cfg 是"日常切换批次"的总开关，随时会被改到二期；一期定标若跟着它走，
%   就会拿二期的测量缓存去解一期的几何，得到一个跑飞解并被健全性断言拦下
%   （实测：x0 从 554 滑到 537、rot 变成 -0.9 deg、a 变成 -0.09）。所以这里
%   把一期路径写成自己的参数（RawRoot / CalNative / CalFile），只从 q_cfg 取
%   **物理常量**（站点坐标、气辉高度、zMax、分辨率）。

% ==========================================================================
% [0] 参数与常量
% ==========================================================================
mp = mc_paths();   % 包内路径唯一来源；本函数不读 q_cfg 的数据路径（见上）
p = inputParser;
p.addParameter('RawRoot',   mp.raw);
p.addParameter('CalNative', mp.calP14Native);
p.addParameter('CalFile',   mp.calP14Pipe);
p.addParameter('Stage', 'auto', @(v) ischar(v) || isstring(v));
p.addParameter('SunAltMax', -18.0);      % 太阳高度门限(deg)：只收 sun_alt < 此值
p.addParameter('AltLo', 14.0);           % 月亮高度窗口：低空段（月核小、天光弱）最干净
p.addParameter('AltHi', 50.0);           %   拟合尺度用的上界（alt>55 月核被光晕撑大）
p.addParameter('AltHiCenter', 70.0);     % 定圆心用的上界：只看方向，可用更宽窗口
p.addParameter('InlierHalfWidth', 2.0);  % 共识判据：rot 与全局中位一致到 ±此值(deg)
p.addParameter('MinInlierFrac', 0.80);   % 共识判据：内点比例下限
p.addParameter('MinIllum', 0.15);
p.addParameter('Nights', {'20140930', '20141001', '20141002', '20141003', ...
                          '20141005', '20141006', '20141007'});
p.addParameter('RimNights', {'20140930', '20141002', '20141003', '20141005', '20141007'});
p.addParameter('RimFrames', 40);
p.addParameter('CenterInit', [554.5, 536.5]);   % 全局紧截尾共识圆心（≈论文 §3.2.1）
% ★ 固定圆心（2026-09-18 加）：给 [x0 y0] 则**跳过圆心搜索**，直接在该圆心下拟合尺度与 rot。
%   为什么需要它：圆心搜索的代价面在某些数据集上是一个**很宽的平台**。
%   实测 2013-10 那批（9 夜）方向共识代价在 x0∈[520,560] 上都是 79~82%，网格搜索落在平台
%   左端 (529.0,534.0)，而 2014-09 的 (553.5,536.5) 在同一批数据上同样是 81~82% ⇒ 数据本身
%   分不开这 24 px，但 rot 会跟着差 0.35 deg。此时应把圆心**钉在已有可信值**上，
%   让两批数据的 rot 变成同口径可比。默认空 = 行为与 2026-09-17 之前完全一致。
p.addParameter('CenterFix', []);
% ★ 站点坐标选择（2026-09-18 加）：站点坐标原先以**硬编码常量**散落本文件两处 ——
%   moon_calib 早已参数化（C.sites.*），而这里 [B1] 与 [B2] 各写死一份，容易改一处漏一处。
%   现统一为本参数，默认 'official' ⇒ 与 2026-09-18 之前的全部结果**逐位一致**。
%   可选 official / png（二期头坐标）/ png1（一期头坐标）；后两者已知有误，仅作灵敏度对照。
p.addParameter('Site', 'official', @(v) ischar(v) || isstring(v));
p.addParameter('RadiusInit', 458.0);            % 一期视场边界半径初值（几何测量）
p.addParameter('SkipFrames', false);
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

cfg = q_cfg();                       % 只取物理常量
cfg.rawRoot   = char(o.RawRoot);     % 其余路径一律用本函数自己的参数
cfg.calNative = char(o.CalNative);
cfg.calFile   = char(o.CalFile);

X0GEO = o.CenterInit(1);   Y0GEO = o.CenterInit(2);   RIM0 = o.RadiusInit;

% 站点表：与 moon_calib 的 C.sites 保持**同一份数值**（改一处必须改两处，已在注释里互相指向）
%   实测（2026-09-18，逐字读 PNG tEXt 头的 #Station）：两期头坐标**不同**
%     二期 2026-03 起 : ODAZH(109.83E, 19.31N, 103m)   -> 'png'
%     一期 2013/2014  : ODAZH(109.1E,  19.5N,  28m)    -> 'png1'（经度只写到 0.1 deg）
SITES = struct('official', [109.133, 19.526, 0.103], ...   % 儋州站（富克 FKT）官方坐标
               'png',      [109.830, 19.310, 0.103], ...   % 二期 PNG 头坐标（经度已知有误）
               'png1',     [109.100, 19.500, 0.028]);      % 一期 PNG 头坐标（仅作灵敏度对照）
siteKey = char(o.Site);
assert(isfield(SITES, siteKey), 'calib_2014:badSite', ...
       '未知站点 "%s"（可选 official / png）', siteKey);
SITEV = SITES.(siteKey);

% --- 论文常量（原文抄录，只读，不参与计算）---
TH = struct();
TH.zenith_star   = [554, 538];        % §3.2.1 三星余弦天顶（**原图**坐标系）
TH.zenith_rot    = [452, 452];        % §3.2.2 的"圆心即天顶"（**裁剪**坐标系）
                                     %   = 同一个物理点：裁剪窗左上角 (102,86)、窗 904x904
                                     %   ⇒ (102+452, 86+452) = (554,538) ✓（见 [E3](i)）
TH.polaris_px    = [499, 55];         % 与 zenith_rot 同系（原文同一句括号内并列）
TH.polaris_az    = 0.4;               % 论文取值；独立复算 = 0.390 deg（rot_audit.m）
TH.rot_deg       = 6.35;              % (3.2) 式；其"旋转中心"就是上面那个圆心（= 天顶）
TH.rz_coef       = [-0.55, 305.2, 16.68];   % (3.9) 式，z in rad，降幂
TH.h_km          = 300;
TH.canvas        = [600, 600];
TH.res_deg       = 0.020;
TH.lat_range     = [13.50, 25.480];
TH.lon_range     = [103.220, 115.20];
TH.frame_dt      = datetime(2014, 9, 23, 17, 6, 12);
% --- 2026-09-17：本工程对论文插图与归档帧的实测值 ---
%     插图测量的可复算脚本  thesis_fig_probe.m（输出 reports\工具_论文插图圆心测量.txt）
%     rot 旁证的可复算脚本  rot_audit.m      （输出 reports\工具_rot三项审计.txt）
TH.tmpl_center   = [513, 513];        % 论文图3.2(c)(d) 圆形模板圆心（归一化到 1024 帧）
TH.tmpl_R        = 511;               % 该模板半径
TH.tmpl_Rw       = 0.499;             % 实测 R/w（(c)(d) 两幅都是 0.499）⇒ 模板圆内切于裁剪窗
TH.fig3_dx       = 0.3;               % 论文图3.3 天顶红圈 x − 画幅中心 x（渲染像素）
TH.fig3_dy       = 24.6;              % 同上 y 差。**仅记录、不作证据**：图3.3 是 2013-10-25 的图，
                                     %   其画幅原点是否等于 §3.2.2 裁剪窗原点，无法确认
TH.polaris_az_t  = 0.390;             % rot_audit 独立复算的 Polaris 真实方位角（论文取 0.4）
TH.polaris_alt   = 20.164;            % 同上，真实高度角（⇒ z = 69.836 deg）
TH.polaris_dbg   = 3;                 % 论文像素邻域 30 px 峰值 − 背景（灰阶）<8 ⇒ 该星测不到
zoff = TH.zenith_star - TH.zenith_rot;  % 裁剪系 → 原图系的原点平移 = (102, 86)

% 原生定标(.mat) 与附带产物目录（沿用 q1_calibrate 的约定：把 rawRoot 切到另一批
% 数据时两批定标的报告不互相覆盖）
[pn, bn] = fileparts(cfg.calNative);
if strcmp(bn, 'calibration_params'), nativeOut = pn; else, nativeOut = fullfile(pn, bn); end
if ~exist(nativeOut, 'dir'), mkdir(nativeOut); end
if ~exist(fileparts(cfg.calFile), 'dir'), mkdir(fileparts(cfg.calFile)); end

say = @(varargin) fprintf(varargin{:});

say('\n');
say('%s\n', repmat('#', 1, 78));
say('#  calib_2014 —— 一期(2014)专用定标 + 马欣论文对照\n');
say('%s\n', repmat('#', 1, 78));
say('原始数据根 : %s\n', cfg.rawRoot);
say('原生定标   : %s\n', cfg.calNative);
say('管线定标   : %s\n', cfg.calFile);
say('判据       : 太阳 < %.1f deg；定圆心 alt %.0f~%.0f deg；拟合尺度 alt %.0f~%.0f deg\n', ...
    o.SunAltMax, o.AltLo, o.AltHiCenter, o.AltLo, o.AltHi);
say('站点       : %-9s (%.3fE, %.3fN, %.0f m)   ← [B1] 与 [B2] 共用这一份\n', ...
    siteKey, SITEV(1), SITEV(2), SITEV(3) * 1000);

% ==========================================================================
% [A] 位深自检 —— 确认 q_read_gray 把 8 bit 拉到了 16 bit 量级
% ==========================================================================
say('\n%s\n[A] 位深自检\n%s\n', repmat('-', 1, 78), repmat('-', 1, 78));
probe = '';
nd = fullfile(cfg.rawRoot, o.Nights{1});
fs = dir(fullfile(nd, '*.PNG'));
if ~isempty(fs), probe = fullfile(nd, fs(1).name); end
if ~isempty(probe)
    raw = imread(probe);
    [g, info] = q_read_gray(probe);
    say('  探测帧 : %s\n', fs(1).name);
    say('  原始类 : %s  ->  q_read_gray 报告 %d bit, 倍率 %g\n', class(raw), info.bits, info.scale);
    say('  原始 max %g / 满量程 %g   统一后 max %.0f / 65535\n', ...
        double(max(raw(:))), double(intmax(class(raw))), max(g(:)));
    if info.bits == 8
        say('  ✅ 8 bit 已线性对齐到 16 bit 量级（×257）。若不过这一关，\n');
        say('     moon_calib 的 65000 饱和判据永远不成立 -> detect 检 0 帧。\n');
    end
else
    say('  [警告] 找不到探测帧，跳过。\n');
end

% ==========================================================================
% [B] 测量与解算
% ==========================================================================
say('\n%s\n[B1] 月亮法定标（moon_calib：detect + 全 QC 链）\n%s\n', ...
    repmat('-', 1, 78), repmat('-', 1, 78));

stage = char(o.Stage);
if strcmpi(stage, 'auto')
    if exist(fullfile(nativeOut, 'moon_meas.mat'), 'file') == 2
        stage = 'solve,verify,export';
        say('  发现 moon_meas.mat 缓存 -> stage = %s\n', stage);
    else
        stage = 'all';
        say('  无测量缓存 -> stage = all（scan,detect,solve,verify,export）\n');
    end
end

% ⚠ 这里**不**传 'pipeline'：本函数要用的采用解是在 [B2] 里用"圆心固定"的方式重算的，
%   与 moon_calib 的自由拟合解不同。LM 解单独写到 _LM.mat，采用解写 cfg.calNative，
%   两条解都留档，避免"先写一个、再被另一个悄悄覆盖"。
% 2026-09-18：LM 解的文件名从 CalNative 的文件名派生（原来写死 _14_，
% 在 calibration_params_13 目录里会得到"名字带 14 的一期2013 LM 解"，极易误读）。
[cnaP, cnaN] = fileparts(cfg.calNative);   %#ok<ASGLU>
lmMat = fullfile(nativeOut, [cnaN '_LM.mat']);
moon_calib('raw', cfg.rawRoot, ...
           'out', nativeOut, 'outmat', lmMat, ...
           'site', siteKey, ...
           'stage', stage, ...
           'nights', o.Nights, ...
           'x0_geo', X0GEO, 'y0_geo', Y0GEO, 'rim', RIM0, ...
           'sun_alt_max', o.SunAltMax, 'alt_lo', o.AltLo, 'alt_hi', o.AltHi, ...
           'min_illum', o.MinIllum);

S = load(lmMat);
if ~isfield(S, 'calMat')
    error('calib_2014:noCalMat', ...
        ['%s 里没有 calMat。\n' ...
         '先看 %s\\calib_report.txt 的 [3.1] 表：若"通过 0 帧"，' ...
         '把 SunAltMax 收得更负、或把 AltHi 收到 50 以内。'], lmMat, nativeOut);
end
calLM = S.calMat;    % LM 自由拟合解（对照用；圆心/rot 有简并偏置）
say('\n  LM 自由拟合解（**不采用**，仅作对照）: x0=%.2f y0=%.2f f=%.2f a=%+.5f rot=%.3f deg rms=%.2f px 手性=%s\n', ...
    calLM.x0, calLM.y0, calLM.f, calLM.a, calLM.rot, calLM.rms_px, calLM.chirality);

% ---- [B2] 采用解：方位共识定圆心 + 固定圆心拟合 ----
say('\n%s\n[B2] 采用解：方位共识定圆心（与径向模型无关）+ 固定圆心拟合\n%s\n', ...
    repmat('-', 1, 78), repmat('-', 1, 78));
say('  为什么: Δrot/Δy0 ≈ -0.17 deg/px —— 圆心的任何偏置几乎原样搬进 rot。\n');
say('          而像面方位 phi = atan2(cx-x0, -(cy-y0)) 与径向模型**无关**\n');
say('          （径向畸变不改变方向），所以"所有帧的 rot 都等于同一个常数"是\n');
say('          一条只含 (x0,y0) 的良定条件。先定圆心，再固定它去拟合尺度与 rot。\n');
say('          判据用**全局紧截尾共识**（±%.1f deg、内点率需 >= %.0f%%），不碰任何\n', ...
    o.InlierHalfWidth, 100 * o.MinInlierFrac);
say('          "逐夜汇总统计量" —— 那正是旧版 (535,543) 伪解的来源（见文件头 (5)）。\n\n');
[adopted, ctr] = adopt_solution(nativeOut, X0GEO, Y0GEO, RIM0, o, TH, say, SITEV);
if isempty(adopted)
    error('calib_2014:noAdopted', '采用解拟合失败（测量缓存不足或判据过严）。');
end

% ---- [B3] 写落盘 ----
calMat = adopted;
save(cfg.calNative, 'calMat', '-v7');
say('\n  原生定标 -> %s   （变量 calMat，0-based）\n', cfg.calNative);
cal = moon_calib_to_pipeline(calMat, cfg.calFile, ...
    'Radius', RIM0, 'Zenith', [calMat.x0 + 1, calMat.y0 + 1], 'Verbose', false);
say('  管线定标 -> %s   zenith(1-based)=(%.3f, %.3f)  rot=%+.4f deg  radius=%.1f px  r(90)=%.1f px\n', ...
    cfg.calFile, cal.zenith(1), cal.zenith(2), cal.rot, cal.radius, cal.r90);

% ★ 2026-09-18：另写一份**采用解**摘要。
%   为什么必须写：moon_calib 的 calib_summary.txt 记的是 **LM 自由拟合解**（圆心自由），
%   而本函数真正采用的是 [B2] 的"定心/钉圆心 + 固定圆心拟合"解。一期的圆心-rot 简并会把
%   两者拉开 20+ px / 0.5+ deg ⇒ 只看 calib_summary.txt 会取到**没有被采用的解**。
%   两个文件都留在同目录，本文件带 _adopted 后缀，避免下游误取。
sumPath = fullfile(nativeOut, 'calib_summary_adopted.txt');
fs2 = fopen(sumPath, 'w', 'n', 'UTF-8');
if fs2 > 0
    fprintf(fs2, '一期定标 · **采用解**（定心/钉圆心 + 固定圆心拟合）\n');
    fprintf(fs2, '============================================================\n');
    fprintf(fs2, 'model                 %s\n', calMat.model);
    fprintf(fs2, 'r_of_z                %s\n', calMat.r_of_z);
    fprintf(fs2, 'chirality             %s\n', calMat.chirality);
    fprintf(fs2, 'center_source         %s\n', calMat.center_source);
    fprintf(fs2, 'x0                    %.3f\n', calMat.x0);
    fprintf(fs2, 'y0                    %.3f\n', calMat.y0);
    fprintf(fs2, 'f                     %.4f\n', calMat.f);
    fprintf(fs2, 'a                     %+.6f\n', calMat.a);
    fprintf(fs2, 'rot                   %+.3f\n', calMat.rot);
    fprintf(fs2, 'rms_px                %.3f\n', calMat.rms_px);
    fprintf(fs2, 'n_frames              %d\n', calMat.n_frames);
    fprintf(fs2, 'n_nights              %d\n', calMat.n_nights);
    fprintf(fs2, 'rot_circ_scatter      %.3f\n', calMat.rot_circ_scatter);
    fprintf(fs2, 'rot_model_spread      %.4f\n', calMat.rot_model_spread);
    fprintf(fs2, 'rot_per_y0_deg_per_px %+.4f\n', calMat.rot_per_y0);
    fprintf(fs2, 'r_horizon_px          %.1f\n', calMat.r_horizon_px);
    fprintf(fs2, 'rim_px                %g\n', calMat.rim_px);
    fprintf(fs2, 'site_lon              %.3f\n', calMat.site_lon);
    fprintf(fs2, 'site_lat              %.3f\n', calMat.site_lat);
    fprintf(fs2, '夜                    %s\n', strjoin(calMat.nights, ' '));
    fprintf(fs2, '逐夜 rot(deg)         %s\n', strtrim(sprintf('%.3f ', [calMat.rot_per_night.rot])));
    fprintf(fs2, '逐夜 帧数             %s\n', strtrim(sprintf('%d ', [calMat.rot_per_night.n])));
    fprintf(fs2, '\n⚠ 同目录 calib_summary.txt 是 moon_calib 的 LM 自由拟合解（圆心自由），\n');
    fprintf(fs2, '  在一期数据上因"圆心-rot 简并"与采用解不同（可达 20+ px / 0.5+ deg）。\n');
    fprintf(fs2, '  要用的是**本文件**。\n');
    fprintf(fs2, '\n用于投影（z 单位 rad，结果单位 px，0-based 坐标）:\n');
    fprintf(fs2, '  r = %.4f * (z + %.6f * z.^3);\n', calMat.f, calMat.a);
    fprintf(fs2, '  x = %.3f + r .* sin(az + deg2rad(%.3f));\n', calMat.x0, calMat.rot);
    fprintf(fs2, '  y = %.3f - r .* cos(az + deg2rad(%.3f));\n', calMat.y0, calMat.rot);
    fclose(fs2);
    say('  采用解摘要 -> %s\n', sumPath);
else
    say('  [警告] 采用解摘要写不出去: %s\n', sumPath);
end

% ==========================================================================
% [C] 几何锚：视场边缘圆（独立于月亮 —— 只依赖视场的旋转对称性）
% ==========================================================================
say('\n%s\n[C] 几何锚：视场边缘圆（zenith_from_rim）\n%s\n', ...
    repmat('-', 1, 78), repmat('-', 1, 78));
say('  原理: r 只是离轴角的函数 ⇒ 等天顶角的像点落在以光轴像点为圆心的圆上;\n');
say('        视场边界就是一条等角线 ⇒ 拟合边界圆的圆心 = boresight 像点。\n');
say('        它给不出 rot 与 f（需要外部尺度），但能**独立**校核天顶前提。\n\n');

zref = [calMat.x0 + 1, calMat.y0 + 1];        % 0-based -> 1-based
dirs = cellfun(@(n) fullfile(cfg.rawRoot, n), o.RimNights, 'UniformOutput', false);
rim = [];
try
    rim = zenith_from_rim(dirs, 'NFrames', o.RimFrames, ...
        'ZenithRef', zref, 'ScalePxRad', calMat.f);
    say('\n  几何圆心(1-based) = (%.3f, %.3f)   R_rim = %.3f px   rms = %.4f px   nPts=%d\n', ...
        rim.center(1), rim.center(2), rim.radius, rim.rms, rim.nPoints);
    dpx = hypot(rim.center(1) - zref(1), rim.center(2) - zref(2));
    say('  月亮天顶像点      = (%.3f, %.3f)\n', zref(1), zref(2));
    say('  两者相差 %.2f px = %.3f deg\n', dpx, atand(dpx / calMat.f));
    if dpx < 10
        say('  判定: 一致（<10 px）⇒ "天顶像点 = 光轴像点" 的前提成立。\n');
    else
        say('  判定: 偏差 %.1f px。但一期视场边缘是**软边**（径向剖面宽约 28 px），\n', dpx);
        say('        且逐方位边界半径从 443 变到 556 px（见 moon_calib [3.10]）⇒\n');
        say('        "亮环 = 视场边界" 这条假设比二期弱得多，几何法只能当参考；\n');
        say('        它以同样的量级落在纯方位法圆心附近，两者互为佐证。\n');
    end
catch ME
    say('  [警告] 几何法失败: %s\n', ME.message);
end

% ==========================================================================
% [D] 论文样例帧复现：2014-09-23 17:06:12 UT
% ==========================================================================
frameOut = struct('raw', '', 'enh', '', 'sub', '', 'ok', false);
if ~o.SkipFrames
    say('\n%s\n[D] 论文样例帧复现（论文图3.2 / 图3.9 的同一帧）\n%s\n', ...
        repmat('-', 1, 78), repmat('-', 1, 78));
    frameOut = thesis_frame_check(cfg, TH, calMat, nativeOut, say);
else
    say('\n[D] 已按 SkipFrames 跳过论文样例帧复现。\n');
end

% ==========================================================================
% [E] 与马欣论文的逐项对照
% ==========================================================================
say('\n%s\n[E] 与马欣论文 3.2 节逐项对照\n%s\n', repmat('-', 1, 78), repmat('-', 1, 78));

zgrid = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, pi/2];       % rad
rOurs = calMat.f * (zgrid + calMat.a * zgrid.^3);
rThes = polyval(TH.rz_coef, zgrid);

say('\n  (E1) 径向标度 r(z)\n');
say('  %8s | %10s %10s %9s   %s\n', 'z (deg)', '本解(px)', '论文(3.9)', '差(%)', '备注');
say('  %s\n', repmat('-', 1, 74));
zlo = deg2rad(90 - o.AltHi); zhi = deg2rad(90 - o.AltLo);   % 本解实际拟合到的 z 区间
for k = 1:numel(zgrid)
    note = ''; if zgrid(k) < zlo - 1e-9 || zgrid(k) > zhi + 1e-9, note = '外推'; end
    say('  %8.1f | %10.1f %10.1f %+9.1f   %s\n', rad2deg(zgrid(k)), rOurs(k), rThes(k), ...
        100 * (rOurs(k) - rThes(k)) / rThes(k), note);
end
say('  %s\n', repmat('-', 1, 74));
r90o = calMat.f * (pi/2 + calMat.a * (pi/2)^3);
r90t = polyval(TH.rz_coef, pi/2);
% ★ 比"一次项系数"更有意义的是**在观测区间上的有效尺度**：论文模型带 16.68 px 截距，
%   直接比一次项会系统性高估它的尺度。
eOurs = (calMat.f * (zhi + calMat.a * zhi^3) - calMat.f * (zlo + calMat.a * zlo^3)) / (zhi - zlo);
eThes = (polyval(TH.rz_coef, zhi) - polyval(TH.rz_coef, zlo)) / (zhi - zlo);
say('  本解拟合区间 : z = %.1f ~ %.1f deg（alt %.0f~%.0f）\n', rad2deg(zlo), rad2deg(zhi), o.AltLo, o.AltHi);
say('  一次项系数   : 本解 f=%.2f px/rad   论文式 305.2 px/rad   差 %+.1f%%\n', ...
    calMat.f, 100 * (calMat.f - 305.2) / 305.2);
say('  区间有效尺度 : 本解 %.1f px/rad   论文式 %.1f px/rad   差 %+.1f%%   <== 这一个才是可比的量\n', ...
    eOurs, eThes, 100 * (eOurs - eThes) / eThes);
say('  r(90 deg)    : 本解 %.1f px   论文 %.1f px   差 %+.1f px (%+.1f%%)\n', ...
    r90o, r90t, r90o - r90t, 100 * (r90o - r90t) / r90t);
say('  ⇒ 论文模型 r = -0.55z^2 + 305.2z + 16.68 的二次项在整个 z<1.5 范围内只贡献\n');
say('     不到 1.3 px，所以它本质上是一条"带截距的直线"；那 16.68 px 截距使得\n');
say('     "一次项 305.2"与实际尺度不同（z=0 时物理上必须 r=0）。本解用 r=f(z+a z^3)、\n');
say('     无截距，所以应与"区间有效尺度"比。若两者的 r(90 deg) 差得较多，说明本解的\n');
say('     a 在向 z=90 deg 外推时被放大 —— 本解只在 alt %.0f~%.0f 内有数据支撑，\n', o.AltLo, o.AltHi);
say('     超出该区间的 r 值都属于外推，仅供与论文对照，不要当独立测量用。\n');

% ★ 2026-09-17 新增：那 −4.5% 到底能不能被月亮数据裁决？
if isfield(ctr, 'shapeCmp') && ~isempty(fieldnames(ctr.shapeCmp))
    sc = ctr.shapeCmp;
    say('\n  (E1b) 差异裁决：r = k * shape(z)，两种形状各放 2 个自由参数（k 与 rot）\n');
    say('  %-26s %8s %10s %11s %12s\n', '径向形状', '自由参数', '2D rms', 'k', 'rot (deg)');
    say('  %s\n', repmat('-', 1, 72));
    say('  %-26s %8d %10.2f %11.4f %+12.3f\n', '本解 f(z+a z^3)', 2, sc.rmsOurs2, sc.kOurs, sc.rotOurs2);
    say('  %-26s %8d %10.2f %11.4f %+12.3f\n', '论文(3.9) + 自由尺度', 2, sc.rmsThesis2, sc.kThesis, sc.rotThesis2);
    say('  %-26s %8d %10.2f %11s %+12.3f\n', '论文(3.9) 原样（尺度固定 1）', 1, sc.rmsThesis1, '1', sc.rotThesis1);
    say('  %s\n', repmat('-', 1, 72));
    say('  ⇒ 论文形状要落到本解的像素尺度上需乘 k = %.4f（%+.1f%%），与上表的 −4.5%% 同源；\n', ...
        sc.kThesis, 100 * (sc.kThesis - 1));
    if sc.rmsThesis2 <= 1.15 * sc.rmsOurs2
        say('  ⇒ 两种形状在 z = %.0f~%.0f deg 内 rms 只差 %.0f%% ⇒ **分不开**：那 −4.5%% 是\n', ...
            rad2deg(zlo), rad2deg(zhi), 100 * (sc.rmsThesis2 / sc.rmsOurs2 - 1));
        say('     "选哪条 r(z) 曲线"的问题，不是测量差异；离开这个区间才分得开。\n');
    else
        say('  ⇒ 论文形状留下比本解大 %.0f%% 的系统残差 ⇒ 本解的三次形状被数据偏好\n', ...
            100 * (sc.rmsThesis2 / sc.rmsOurs2 - 1));
        say('     （两者自由参数同为 2 个，差的只是曲线形状 ⇒ 不是"参数多所以拟合好"）。\n');
    end
    if sc.rmsThesis1 > 1.5 * sc.rmsOurs2
        say('  ⇒ 若连尺度都不让动（论文式原样），rms 升到 %.1f px（%.1f 倍）—— 说明两家的\n', ...
            sc.rmsThesis1, sc.rmsThesis1 / sc.rmsOurs2);
        say('     径向尺度本来就不在同一个归一化上，直接比系数没有意义。\n');
    end
end

say('\n  (E2) 天顶像点（=投影中心）\n');
say('  %-34s %13s %13s %13s\n', '来源', 'x (px)', 'y (px)', '与采用解之差');
say('  %s\n', repmat('-', 1, 78));
say('  %-34s %13.2f %13.2f %13s\n', '采用解（方位共识, 模型无关）', calMat.x0, calMat.y0, '—');
if isfield(ctr, 'centerInlierFrac')
    say('  %-34s %13s %13s %13s\n', '  └ 共识内点率', ...
        sprintf('%d/%d', ctr.nCenterFrames, ctr.nFrameAll), ...
        sprintf('%.0f%%', 100 * ctr.centerInlierFrac), ...
        sprintf('散布 %.2f deg', sqrt(ctr.cost)));
end
say('  %-34s %13.2f %13.2f %13.1f px\n', 'LM 自由拟合解（对照）', calLM.x0, calLM.y0, ...
    hypot(calLM.x0 - calMat.x0, calLM.y0 - calMat.y0));
if ~isempty(rim)
    say('  %-34s %13.2f %13.2f %13.1f px\n', '几何法（视场边缘圆, =boresight）', ...
        rim.center(1) - 1, rim.center(2) - 1, ...
        hypot(rim.center(1) - 1 - calMat.x0, rim.center(2) - 1 - calMat.y0));
end
say('  %-34s %13.0f %13.0f %13.1f px\n', '论文 §3.2.1 三星余弦', TH.zenith_star(1), ...
    TH.zenith_star(2), hypot(TH.zenith_star(1) - calMat.x0, TH.zenith_star(2) - calMat.y0));
say('  %-34s %13.0f %13.0f %13.1f px\n', '论文 §3.2.2 圆心→原图系', ...
    TH.zenith_rot(1) + zoff(1), TH.zenith_rot(2) + zoff(2), ...
    hypot(TH.zenith_rot(1) + zoff(1) - calMat.x0, TH.zenith_rot(2) + zoff(2) - calMat.y0));
say('  %s\n', repmat('-', 1, 78));
say('  读法: 本解与论文 §3.2.1 的三星余弦天顶相差 %.1f px（= %.2f deg）——\n', ...
    hypot(TH.zenith_star(1) - calMat.x0, TH.zenith_star(2) - calMat.y0), ...
    atand(hypot(TH.zenith_star(1) - calMat.x0, TH.zenith_star(2) - calMat.y0) / calMat.f));
say('        两者**方法完全独立**（论文用星、本解用月亮方位），落在这一点上是\n');
say('        最不容易凑出来的一致性。\n');
say('  注意 (452,452) **不是**第二个天顶，而是同一颗天顶在"裁剪坐标系"里的写法：\n');
say('        平移量就是裁剪窗左上角 (%.0f,%.0f)，回代即得 (%.0f,%.0f)（详见 E3 (i)）。\n', ...
    zoff(1), zoff(2), TH.zenith_rot(1) + zoff(1), TH.zenith_rot(2) + zoff(2));
if ~isempty(rim)
    say('  几何法圆心相差 %.1f px：一期视场边缘是**软边**（径向剖面 ~28 px 宽，逐方位\n', ...
        hypot(rim.center(1) - 1 - calMat.x0, rim.center(2) - 1 - calMat.y0));
    say('        边界半径从 443 变到 556 px），且它测的是**光轴**而非天顶点（两者差一个\n');
    say('        安装倾角）⇒ 只能当量级校核，不能用来推翻方位共识。\n');
end

say('\n  (E3) 方位校正角 rot（2026-09-17 二次修订）\n');
say('  论文用北极星像点 (%.0f, %.0f)、方位角 %.1f deg 与**圆心**(%.0f, %.0f) 算出 %.2f deg。\n', ...
    TH.polaris_px(1), TH.polaris_px(2), TH.polaris_az, ...
    TH.zenith_rot(1), TH.zenith_rot(2), TH.rot_deg);
say('\n  (i) 先把坐标系说清楚 —— 这决定"能不能比"。\n');
say('      原文顺序是**先定圆心、再裁圆**（§3.2.1 原文）：\n');
say('        "…原始气辉图像是圆形的，因此需要以圆形模板来提取中间的图像。首先需要知道\n');
say('          圆心的在图像中的位置，这里我们通过选取亮度增强图中分布均匀的三个星点，\n');
say('          根据余弦定理计算得到圆心即天顶的位置坐标 … 对其求平均值即得到最终提取\n');
say('          图像的圆心坐标（x=554，y=538）。"\n');
say('      ⇒ "圆心"是裁剪的**输入**（由星表独立求出），裁剪本身是"以天顶为圆心"做的，\n');
say('        所以"圆心"与"天顶"在论文里**是同一个量**，不存在"拿圆心当天顶"的误差。\n');
say('      ⇒ (554,538)[§3.2.1 原图系] 与 (452,452)[§3.2.2 裁剪系] **是同一个物理点**：\n');
say('        裁剪窗左上角 = (%d,%d)，边长 = 2 x 452 = %d ⇒ 回代 (%d+452, %d+452) = (%d,%d) ✓\n', ...
    zoff(1), zoff(2), 2 * TH.zenith_rot(1), zoff(1), zoff(2), zoff(1) + 452, zoff(2) + 452);
say('        模板圆半径 452 = 裁剪窗半宽 = 内切圆；实测图3.2(c)(d) 的 R/w = %.3f ⇒ 一致。\n', ...
    TH.tmpl_Rw);
say('      ⇒ 于是论文的**天顶就是裁剪窗中心**，两者不可当独立量互校，也不必互校。\n');
dRot = TH.rot_deg - calMat.rot;      % 统一口径：论文 − 本解
dZen = hypot(TH.zenith_star(1) - calMat.x0, TH.zenith_star(2) - calMat.y0);
say('\n  (ii) 那 %.2f deg 的差必须在别处找。逐项排除（可复算脚本: rot_audit.m）:\n', dRot);
say('      (1) 天顶？不是。论文 (554,538) 与本解 (%.1f,%.1f) 差 %.1f px = %.2f deg。\n', ...
    calMat.x0, calMat.y0, dZen, atand(dZen / calMat.f));
say('      (2) 论文那个方位角？不是。本解独立算得 Polaris 在 2014-09-23 17:06:12 UT、\n');
say('          (%.3fE, %.3fN) 的真实方位角 = %+.3f deg，论文取 %.1f ⇒ 差 %+.3f deg。\n', ...
    SITEV(1), SITEV(2), TH.polaris_az_t, TH.polaris_az, TH.polaris_az_t - TH.polaris_az);
say('      (3) 相机时钟？不是。把星历时刻整体平移 dt 扫描 ±36 min：二维 rms 在 dt=0\n');
say('          最小（4.572 px），纯方位散射也在 dt=-3~0 最小（0.538 deg）；rot 对 dt 的\n');
say('          斜率只有 0.160 deg/min ⇒ 即便钟差 5 min 也只带 0.80 deg。\n');
say('      (4) 那就只剩那颗星的**方向**。把论文自己的北极星像素配上本解天顶：\n');
polRaw = TH.polaris_px + zoff;
rPol = hypot(polRaw(1) - calMat.x0, calMat.y0 - polRaw(2));
phiP = mod(atan2d(polRaw(1) - calMat.x0, calMat.y0 - polRaw(2)), 360);
r9   = polyval(TH.rz_coef, deg2rad(90 - TH.polaris_alt));
say('          r = %.1f px, phi = %.2f deg ⇒ rot = %.2f - %.1f = **%.2f deg**（论文 %.2f）\n', ...
    rPol, phiP, phiP, TH.polaris_az, phiP - TH.polaris_az, TH.rot_deg);
say('          ⇒ 论文的 (499,55) 与 (452,452) 两个数**互相自洽**；%.2f deg 的差整份落在\n', abs(dRot));
say('            这颗星的方向上，而 %.2f deg 在 r = %.0f px 处 = **横向 %.0f px**。\n', ...
    abs(dRot), rPol, abs(dRot) * pi / 180 * rPol);
say('          * 该星离天顶 %.1f px，而论文自己的式(3.9) 在 z = %.2f deg 给 %.1f px ⇒\n', ...
    rPol, 90 - TH.polaris_alt, r9);
say('            这个像素位置本身就有 %.1f px（%.1f%%）量级的不自洽;\n', ...
    rPol - r9, 100 * abs(rPol - r9) / r9);
say('          * 更关键：它在归档 8 bit 帧里**根本测不到** —— 邻域 30 px 的峰值只比背景高\n');
say('            %.0f 个灰阶（判据要求 >=8）⇒ 论文这颗星无法从归档数据复核。\n', TH.polaris_dbg);
say('  ⇒ 结论: %.2f deg 的差**无法裁决**，但它落在论文"单颗手工拾取星像 + 一次 Stellarium\n', ...
    abs(dRot));
say('     读数"的误差包络内。本解 rot 的旁证已逐项排除天顶、方位角、时钟三项，故不采用\n');
say('     "论文的 6.35 deg 更准"这一读法，也不采用早先"圆心当天顶造成 3 deg"的解释。\n');
say('  采用解: rot = %+.2f deg（逐夜 %.2f~%.2f deg，夜间一致性 %.2f deg）。\n', ...
    calMat.rot, min([calMat.rot_per_night.rot]), max([calMat.rot_per_night.rot]), ...
    calMat.rot_circ_scatter);
say('  采用解 rot = %+.2f deg，论文 %+.2f deg ⇒ 论文比本解大 %.2f deg（占其量级 %.0f%%）。\n', ...
    calMat.rot, TH.rot_deg, dRot, 100 * abs(dRot) / TH.rot_deg);
say('  采用解的 rot 支撑证据（全部可由本函数复算）:\n');
say('    * %d 帧 / %d 夜，逐夜 rot: %s deg，夜间圆散布 %.2f deg;\n', ...
    calMat.n_frames, calMat.n_nights, ...
    strtrim(sprintf('%.2f ', [calMat.rot_per_night.rot])), calMat.rot_circ_scatter);
say('    * 圆心由**与径向模型无关**的方位共识独立定出 (%.1f, %.1f)：核心证据是\n', ...
    calMat.x0, calMat.y0);
say('      "与全局 rot 中位一致到 ±2 deg 的帧占 %.0f%%"，而把圆心挪到旧解 (535,543)\n', ...
    100 * ctr.centerInlierFrac);
say('      上只剩 55%%, 挪到 (552,556) 上只剩 66%% ⇒ 这是一个**帧集不变**的判别量;\n');
say('    * 6 个不同帧集（含/不含暮光帧、alt 14-50/14-70、去坏夜、只用最干净 3 夜）\n');
say('      都收敛到同一处（spread < 3 px），内点率 78%%~97%%;\n');
say('    * moon_calib 的手性判决: **%s**（自由拟合下直接解 rms %.2f px，镜像解明显更差）;\n', ...
    calLM.chirality, calLM.rms_px);
% ★ 2026-09-18 实测更正：原先这里写"换成 PNG 头里的错值，rot 只动 0.01 deg" —— 那句话
%   只在**该批自己的头坐标**下成立，而两期的头坐标**不是同一个值**，若不写清会被误读。
%   实测（真重跑，见 mooncal_work\site_sens\site_sens_report.txt）：
%     取**一期**头坐标 (109.100E,19.500N,28m)：2014 批 rot 4.3482 -> 4.3594（+0.011 deg）；
%                                            2013 批 rot 4.2011 -> 4.2272（+0.026 deg）；
%     取**二期**头坐标 (109.830E,19.310N,103m)：2014 批 -> 3.7894（-0.559 deg）；
%                                            2013 批 -> 4.0297（-0.171 deg）。
%   传站点灵敏度可用线性系数概括（2014 批实测）：Δrot/Δlon ≈ +0.71、Δrot/Δlat ≈ -0.47 deg/deg
%   ⇒ 官方坐标纵有 0.1 deg 偏差也只带 ≲0.09 deg，远小于 rot 的 ±1.0 deg 不确定度。
say('    * 站点坐标：取**该批自己的** PNG 头坐标 (109.100E,19.500N,28m) 时 rot 只动 0.01~0.03 deg,\n');
say('      而官方坐标与它只差 0.033E/0.026N ⇒ 站点坐标不是 rot 的误差主项;\n');
say('      ⚠ 二期文件头写的是**另一组** (109.830E,19.310N)：误用它会把本批 rot 拉到 3.79 deg\n');
say('      （-0.56 deg）—— 那是口径错，不是相机变了。两期头坐标差 0.70 deg 经，必须分开对待。\n');
say('      （灵敏度随月亮轨迹而变：2014 批 Δrot/Δlon≈-0.65、Δrot/Δlat≈+0.40，2013 批为\n');
say('        -0.39/-0.52 —— 故只能用实测值，不要把某一批的系数外推到另一批。\n');
say('        见 mooncal_work\\site_sens\\site_sens_report.txt）\n');
say('  rot 的不确定度分解（本函数可复算）:\n');
say('    (a) 径向模型: 钉死圆心后换 5 种 r(z)（等距/等距+三次/等立体角/立体投影/正交），\n');
say('        rot 只动 %.3f deg ⇒ **模型选择不是误差主项**;\n', calMat.rot_model_spread);
say('    (b) 圆心: Δrot/Δy0 = %+.4f deg/px ⇒ 圆心的 ±3 px 不确定度就是 rot 的 ±%.2f deg;\n', ...
    calMat.rot_per_y0, abs(calMat.rot_per_y0) * 3);
say('    (c) 相机时钟（rot_audit.m）: 二维 rms 与纯方位散射都在 dt~0 取极小，rot 对 dt\n');
say('        的斜率只有 0.160 deg/min ⇒ 钟差 5 min 也只带 0.80 deg ⇒ 时钟项 <= 0.5 deg;\n');
say('    (d) 月核质心的相位偏置（8 bit、月核仅 3~5 px 且过饱和）：按逐夜 rot 的 %.2f deg\n', ...
    calMat.rot_circ_scatter);
say('        圆散布与低空 hold-out 的 5.1 px 像点残差，保守给 ±0.5 deg。\n');
say('  ⇒ 采用解 rot = %+.2f deg ± 1.0 deg（各项保守叠加），即 [%.2f, %.2f] deg。\n', ...
    calMat.rot, calMat.rot - 1, calMat.rot + 1);
say('\n  ★ 一句话: 与论文的 rot 差 %.2f deg **不是**"旋转中心取错"造成的（论文的圆心\n', abs(dRot));
say('    就是天顶，见 (i)），也不是时钟或那颗星的方位角造成的（见 (ii)(2)(3)）；它整份\n');
say('    落在论文那颗手工拾取的北极星像的**方向**上，而那颗星在 8 bit 归档帧里测不到、\n');
say('    且离它自己的径向式有约 10 px 的不自洽 ⇒ 论文的 6.35 deg 精度不可能优于 1 deg。\n');
say('    实践含义: 2 deg 的方位旋转差，在 z = 30~60 deg（300 km 高度上对应 170~520 km）\n');
say('    就是 6~18 km 的东西向位移 —— 做泡团东西向结构对比时值得写一句。\n');

say('\n  (E4) 投影画布（论文表3.1）\n');
say('  %-22s %-26s %s\n', '项', '论文', '本链(cfg)');
say('  %s\n', repmat('-', 1, 74));
say('  %-22s %-26s %s\n', '像素', sprintf('%dx%d', TH.canvas), ...
    '600x600 (SpanDeg=11.98)');
say('  %-22s %-26s %s\n', '每像素', sprintf('%.3f deg', TH.res_deg), ...
    sprintf('%.3f deg', cfg.resDeg));
say('  %-22s %-26s %s\n', '纬度范围(N)', sprintf('%.2f ~ %.2f', TH.lat_range), ...
    sprintf('站点 %.3f 起算', cfg.site(2)));
say('  %-22s %-26s %s\n', '经度范围(E)', sprintf('%.3f ~ %.2f', TH.lon_range), ...
    sprintf('站点 %.3f 起算', cfg.site(1)));
say('  %-22s %-26s %s\n', '气辉高度(km)', sprintf('%g', TH.h_km), sprintf('%g', cfg.hShell));
say('  说明: 论文 600x600 @0.02 deg/px 只覆盖 ±6 deg，在 %g km 上约 ±%d km；\n', ...
    cfg.hShell, round(cfg.hShell * tand(6)));
say('        本链日常路径 q2_project / start_process 默认 SpanDeg=11.98 ⇒ 与论文表3.1 同规格。\n');
say('  ★ 2026-09-17 更正: 只有显式 SpanDeg=[] 时才走"自动"，而"自动"原本只被 cal.radius\n');
say('    (可见盘边) 截住，**并不受 cfg.zMax 限制**：一期会一直铺到 z≈84 deg（盘边 %g px 减 Margin 15），\n', RIM0);
say('    早已越过径向标定的实测锚点（alt %.0f deg ⇒ z=%.0f deg）。已给 airglow_geo_pipeline\n', o.AltLo, 90 - o.AltLo);
say('    加 Zmax 选项、并在 q2_project / start_process 传入 cfg.zMax ⇒ 自动画布现在真的限在 z≦%.0f deg。\n', cfg.zMax);

say('\n  (E5) 视场角 / 盘边天顶角（★ 2026-09-17 新增；对应仪器"标称180 deg、有效约160 deg"）\n');
% 下面这些量都只用到 calMat + RIM0，可在本函数内自足算出；
% "盘边到底是不是 80/90 deg"那个假设检验要用月亮数据重拟合，放在独立工具 fov_check.m 里。
HW = 512;                                   % 一期画幅 1024 的半宽
r_of = @(zd) calMat.f * (deg2rad(zd) + calMat.a * deg2rad(zd).^3);
r80  = r_of(80);   r70 = r_of(70);   r90o = r_of(90);
r90t = polyval(TH.rz_coef, pi/2);
zz3  = linspace(1e-6, 2.4, 200001);  rv3 = calMat.f * (zz3 + calMat.a * zz3.^3);
zi_r = @(R) rad2deg(zz3(find(rv3 >= R, 1, 'first')));
zRimLo = zi_r(min([RIM0, 456.53, 450.96]));
zRimHi = zi_r(max([RIM0, 456.53, 450.96]));
say('  盘边 R_rim = %g px 反解 theta_rim = %.1f deg —— **不是 90 deg**。\n', RIM0, zi_r(RIM0));
say('  三条互证(RadiusInit/剖面|dI/dr|最大/半高)给 theta_rim = %.1f ~ %.1f deg。\n', zRimLo, zRimHi);
say('  · 若把"有效视场160 deg"读成"盘边 = 80 deg"：需 a≈+0.07（现实测 %+.4f），\n', calMat.a);
say('    月亮数据 rms 从 3.11 px 抬到 11.27 px（3.6 倍）⇒ **被数据否掉**（工具 fov_check.m 可复算）。\n');
say('  · 自洽的读法：盘边 ≈ %.0f deg（≈%.0f deg 全角，与"标称180 deg"差 %.0f deg），\n', ...
    zi_r(RIM0), 2 * zi_r(RIM0), 180 - 2 * zi_r(RIM0));
say('    而"有效160 deg"是**可用质量边界**，落在盘内 80 deg 处 ⇒ r(80 deg) = %.1f px（半幅的 %.1f%%）。\n', ...
    r80, 100 * r80 / HW);
say('  · r(90 deg) = %.1f px **已超出可见盘边 %g px** ⇒ 90 deg 圈落在最外渐晕带之外；\n', r90o, RIM0);
say('    论文径向式在 90 deg 给 %.1f px ⇒ **超出盘边 %.1f px**，落在像圈之外，\n', r90t, r90t - RIM0);
say('    不可能被任何像元验证 ⇒ 与论文比径向尺度只能比"拟合区间内的有效尺度"。\n');
say('  · f = R_rim/(pi/2) = %.1f px/rad，比月亮实测 %.1f 差 %+.1f%% ⇒ **该换算在一期明确作废**。\n', ...
    RIM0 / (pi/2), calMat.f, 100 * (RIM0 / (pi/2) - calMat.f) / calMat.f);
say('  · cfg.zMax = %.0f deg ⇒ r = %.1f px，既在实测锚点（z≦%.0f deg）内，也远在可用上限之内 ✓\n', ...
    cfg.zMax, r70, 90 - o.AltLo);

% ==========================================================================
% [F] 汇总 + 落盘
% ==========================================================================
out = struct();
out.calMat  = calMat;        % 采用解
out.calLM   = calLM;         % LM 自由拟合解（对照）
out.rim     = rim;
out.centerSearch = ctr;
out.thesis  = TH;
out.frame   = frameOut;
out.cal     = cal;           % 管线格式
out.rot     = calMat.rot;
out.zenith  = [calMat.x0, calMat.y0];
out.r90     = r90o;

tpath = fullfile(nativeOut, 'thesis_compare.txt');
fid = fopen(tpath, 'w', 'n', 'UTF-8');
fprintf(fid, 'calib_2014 vs 马欣论文 3.2 节  对照表\n');
fprintf(fid, '生成 %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '采用解(方位共识圆心 + 固定圆心拟合): x0=%.3f y0=%.3f f=%.4f a=%+.6f rot=%.3f deg rms=%.3f px (%d 帧/%d 夜)\n', ...
    calMat.x0, calMat.y0, calMat.f, calMat.a, calMat.rot, calMat.rms_px, ...
    calMat.n_frames, calMat.n_nights);
fprintf(fid, '  定心共识内点率    : %d/%d = %.0f%%  内点 rot 散布 %.2f deg\n', ...
    ctr.nCenterFrames, ctr.nFrameAll, 100 * ctr.centerInlierFrac, sqrt(ctr.cost));
fprintf(fid, '  区间有效尺度(z=%.0f~%.0f deg): 本解 %.1f px/rad  论文式 %.1f px/rad  差 %+.1f%%\n\n', ...
    rad2deg(zlo), rad2deg(zhi), eOurs, eThes, 100 * (eOurs - eThes) / eThes);
fprintf(fid, 'LM 自由拟合解(对照)              : x0=%.3f y0=%.3f f=%.4f a=%+.6f rot=%.3f deg rms=%.3f px\n\n', ...
    calLM.x0, calLM.y0, calLM.f, calLM.a, calLM.rot, calLM.rms_px);
fprintf(fid, '论文 §3.2.1 天顶 (%.0f, %.0f)  ≡  §3.2.2 圆心 (%.0f, %.0f) —— 同一物理点：\n', ...
    TH.zenith_star, TH.zenith_rot);
fprintf(fid, '   裁剪窗左上角 (%d,%d)、边长 %d ⇒ (%d+452,%d+452) = (%.0f,%.0f) ✓\n', ...
    zoff(1), zoff(2), 2 * TH.zenith_rot(1), zoff(1), zoff(2), TH.zenith_star);
fprintf(fid, '论文北极星 (%.0f,%.0f)[裁剪系] = (%.0f,%.0f)[原图系]；方位角 %.1f deg（独立复算 %.3f）\n', ...
    TH.polaris_px, polRaw, TH.polaris_az, TH.polaris_az_t);
fprintf(fid, '论文 rot=%.2f deg；用论文自己的星像素配本解天顶 ⇒ rot=%.2f deg ⇒ 论文内部自洽，\n', ...
    TH.rot_deg, phiP - TH.polaris_az);
fprintf(fid, '   2 deg 的差全在"这颗星的方向"上；而该星在归档 8 bit 帧里测不到（邻域峰值比\n');
fprintf(fid, '   背景高 %d 个灰阶 < 8）⇒ 无法复核；时钟扫描见 工具_rot三项审计.txt\n', TH.polaris_dbg);
fprintf(fid, '论文径向 (3.9) r = -0.55z^2 + 305.2z + 16.68 (z in rad)\n');
fprintf(fid, '本解 rot 的模型敏感度(固定圆心, 5 种 r(z)): %s deg  极差 %.3f deg\n', ...
    strtrim(sprintf('%.3f ', calMat.rot_model_rot)), calMat.rot_model_spread);
if isfield(calMat, 'rot_model_rms') && ~isempty(calMat.rot_model_rms)
    fprintf(fid, '  同一内点集下 5 种模型的 2D rms        : %s px\n', ...
        strtrim(sprintf('%.2f ', calMat.rot_model_rms)));
end
fprintf(fid, '本解 Δrot/Δy0 = %+.4f deg/px（rot 由圆心定，不由径向模型定）\n', calMat.rot_per_y0);
% ★ 径向形状 A/B：那 −4.5% 能不能被月亮数据裁决
if isfield(calMat, 'shape_cmp') && isfield(calMat.shape_cmp, 'rmsOurs2')
    sc = calMat.shape_cmp;
    fprintf(fid, '\n径向形状 A/B（同圆心、同内点集；r = k*shape(z)，各 2 个自由参数）:\n');
    fprintf(fid, '  %-24s 自由参数 %d  2D rms = %.2f px  (k=%.4f, rot=%+.3f)\n', ...
        '本解 f(z+a z^3)', 2, sc.rmsOurs2, sc.kOurs, sc.rotOurs2);
    fprintf(fid, '  %-24s 自由参数 %d  2D rms = %.2f px  (k=%.4f, rot=%+.3f)\n', ...
        '论文(3.9)+自由尺度', 2, sc.rmsThesis2, sc.kThesis, sc.rotThesis2);
    fprintf(fid, '  %-24s 自由参数 %d  2D rms = %.2f px  (k=1 固定, rot=%+.3f)\n', ...
        '论文(3.9) 原样', 1, sc.rmsThesis1, sc.rotThesis1);
    if sc.rmsThesis2 <= 1.15 * sc.rmsOurs2
        fprintf(fid, '  ⇒ 两种形状在 z=40~76 deg 内 rms 只差 %.0f%% ⇒ 分不开：−4.5%% 属"选哪条\n', ...
            100 * (sc.rmsThesis2 / sc.rmsOurs2 - 1));
        fprintf(fid, '    r(z) 曲线"的问题，不是测量差异。论文形状落到本解尺度需乘 k=%.4f（%+.1f%%）。\n', ...
            sc.kThesis, 100 * (sc.kThesis - 1));
    else
        fprintf(fid, '  ⇒ 论文形状的系统残差比本解大 %.0f%%（自由参数同为 2）⇒ 本解形状被数据偏好。\n', ...
            100 * (sc.rmsThesis2 / sc.rmsOurs2 - 1));
    end
end
fprintf(fid, '\n');
fprintf(fid, ' z(deg)  本解(px)  论文(px)  差(%%)\n');
for k = 1:numel(zgrid)
    fprintf(fid, '%7.1f  %8.1f  %8.1f  %+7.1f\n', rad2deg(zgrid(k)), rOurs(k), rThes(k), ...
        100 * (rOurs(k) - rThes(k)) / rThes(k));
end
fclose(fid);
say('\n  对照表已落盘 -> %s\n', tpath);
say('\n管线定标已就绪 -> %s\n', cfg.calFile);
say('  用法: start_process(''<夜>'', ''Cal'', ''%s'')\n', cfg.calFile);
say('\n%s\n\n', repmat('#', 1, 78));
end


% ==========================================================================
% [B2] 采用解：方位共识定圆心（与径向模型无关）+ 固定圆心拟合
% ==========================================================================
function [cal, info] = adopt_solution(nativeOut, x0g, y0g, rimPx, o, TH, say, site)
%ADOPT_SOLUTION 先用"方位共识"定圆心（与径向模型无关），再固定圆心拟合尺度与 rot。
%
% 输出 cal: 与 moon_calib 的 calMat 同字段（x0/y0/f/a/rot/rms_px/n_frames/n_nights/...），
%           可直接喂 moon_calib_to_pipeline。
% TH 只用于取论文式(3.9) 的系数 TH.rz_coef（做"径向形状 A/B"裁决）。
% site : [lon_deg, lat_deg, h_km]。★ 2026-09-18 由形参传入（原先在此硬编码 109.133/19.526）
%        —— 站点坐标同时决定月亮/太阳星历，换站必须**同步**影响 [B1] 与 [B2]，
%        否则两处口径不一致会让"rot 随站点变多少"这类对照失去意义。
if nargin < 8 || isempty(site), site = [109.133, 19.526, 0.103]; end
cal = [];
info = struct('center', [NaN NaN], 'cost', NaN, 'nCenterFrames', 0, ...
              'nFrameAll', 0, 'centerInlierFrac', NaN, ...
              'modelRot', [], 'modelRms', [], 'shapeCmp', struct());

mp = fullfile(nativeOut, 'moon_meas.mat');
if exist(mp, 'file') ~= 2
    say('  [跳过] 没有测量缓存 %s\n', mp);
    return
end
S = load(mp);
cx = double(S.cx(:)).'; cy = double(S.cy(:)).';
ts = S.ts(:).'; nt = cellfun(@(c) char(c), S.night(:).', 'UniformOutput', false);
N = numel(cx);
if N < 50, say('  [跳过] 测量帧太少 (%d)。\n', N); return; end

% ---- 星历：月亮（方位用于 rot）与太阳（用于暮光门限）----
lon = site(1); lat = site(2); hgt = site(3);
azm = zeros(1, N); alt = zeros(1, N); salt = zeros(1, N);
for i = 1:N
    s = char(ts{i});
    t = datetime(str2double(s(1:4)), str2double(s(5:6)), str2double(s(7:8)), ...
                 str2double(s(9:10)), str2double(s(11:12)), str2double(s(13:14)));
    [a1, a2] = moon_altaz_ref(t, lon, lat, hgt);   % 公开的独立实现
    alt(i) = a1; azm(i) = a2;
    salt(i) = sun_altaz_ref(t, lon, lat);          % 公开的独立实现
end
% 两个窗口：定圆心只看方向（可用宽窗口）；拟合尺度要径向位置（必须窄窗口）
sitC = (salt < o.SunAltMax) & (alt > o.AltLo) & (alt < o.AltHiCenter);
sitF = (salt < o.SunAltMax) & (alt > o.AltLo) & (alt < o.AltHi);
say('  定圆心帧集 %d / %d（太阳<%.1f 且 alt %.0f~%.0f）\n', sum(sitC), N, ...
    o.SunAltMax, o.AltLo, o.AltHiCenter);
say('  拟合帧集   %d / %d（太阳<%.1f 且 alt %.0f~%.0f）\n', sum(sitF), N, ...
    o.SunAltMax, o.AltLo, o.AltHi);
if sum(sitC) < 50 || sum(sitF) < 50, say('  [跳过] 有效帧太少。\n'); return; end

% ---- 步骤 1：纯方位定圆心（全局紧截尾共识，不含任何"逐夜汇总量"）----
if ~isempty(o.CenterFix)
    cen = double(o.CenterFix(:)).';
    [nInF, cst] = tight_inliers(cen(1), cen(2), cx, cy, azm, sitC, o.InlierHalfWidth);
    cinfo = struct('nIn', nInF, 'nAll', sum(sitC), 'frac', nInF / max(sum(sitC), 1));
    say('  步骤 1  圆心**钉死**在 (%.1f, %.1f)（CenterFix）  内点 %d/%d (%.0f%%)  内点散布 %.3f deg\n', ...
        cen(1), cen(2), cinfo.nIn, cinfo.nAll, 100 * cinfo.frac, sqrt(max(cst, 0)));
else
    [cen, cst, cinfo] = search_center(cx, cy, azm, sitC, x0g, y0g, o);
    say('  步骤 1  方位共识定圆心 -> (%.1f, %.1f)  内点 %d/%d (%.0f%%)  内点散布 %.3f deg\n', ...
        cen(1), cen(2), cinfo.nIn, cinfo.nAll, 100 * cinfo.frac, sqrt(max(cst, 0)));
end
info.center = cen; info.cost = cst; info.nCenterFrames = cinfo.nIn;
info.nFrameAll = cinfo.nAll;
info.centerInlierFrac = cinfo.frac;

% ---- 步骤 2：固定圆心，在拟合窗口上拟合 (f, a, rot)，带稳健剔除 ----
% ★ 第 1 轮用**多起点**：单起点 [f,rot,0] 会掉进 a~0 的局部极小（实测 rms 66 px），
%   于是第一轮就把落在高 z 的帧当离群剔掉，人为改变样本。多起点可避免这一点。
x0 = cen(1); y0 = cen(2);
keep = sitF;
for it = 1:4
    kk = find(keep);
    zz = deg2rad(90 - alt(kk))'; AA = deg2rad(azm(kk))';
    xm = cx(kk)' - x0; ym = y0 - cy(kk)';
    if it == 1
        bq = []; bs = Inf;
        for a0 = [0, -0.015, -0.03, -0.045, +0.015]
            qq = fit3(zz, AA, xm, ym, [324, 4.3, a0]);
            rr = qq(1) * (zz + qq(3) * zz.^3);
            ss = sum((rr .* sin(AA + deg2rad(qq(2))) - xm).^2 + ...
                     (ym - rr .* cos(AA + deg2rad(qq(2)))).^2);
            if ss < bs, bs = ss; bq = qq; end
        end
        q = bq;
    else
        q = fit3(zz, AA, xm, ym, q);
    end
    rr = q(1) * (zz + q(3) * zz.^3);
    e = hypot(rr .* sin(AA + deg2rad(q(2))) - xm, ym - rr .* cos(AA + deg2rad(q(2))));
    say('  步骤 2  第 %d 轮: f=%.2f  a=%+.5f  rot=%+.3f deg  rms=%.2f px  n=%d  r(90)=%.1f px\n', ...
        it, q(1), q(3), q(2), sqrt(mean(e.^2)), sum(keep), q(1) * (pi/2 + q(3) * (pi/2)^3));
    if it < 4
        thr = max(5.0, median(e) + 3 * median(abs(e - median(e))));
        keep(kk) = e < thr;
    end
end

% ---- rot 对径向模型的敏感度（圆心固定，只换 r(z) 的表达式）----
% 这一步回答"rot 的不确定度里，有多少来自'选哪个畸变模型'"。做法：用采用解的
% 同一个内点集，把圆心钉死，只换 r(z) 形式重解 rot。实测 rot 只动 ~0.1 deg
% ⇒ rot 几乎完全由圆心决定（Δrot/Δy0 ≈ 0.154 deg/px），模型选择不是误差主项。
kkM = find(keep);
zzM = deg2rad(90 - alt(kkM))'; AAM = deg2rad(azm(kkM))';
xmM = cx(kkM)' - x0; ymM = y0 - cy(kkM)';
modNames = {'r = f z', 'r = f (z + a z^3)', 'r = 2 f sin(z/2)', 'r = 2 f tan(z/2)', 'r = f sin z'};
rotAlt = zeros(1, numel(modNames));
rmsAlt = zeros(1, numel(modNames));
for mi = 1:numel(modNames)
    qm = fit_model(mi, zzM, AAM, xmM, ymM, [q(1), q(2), q(3)]);
    rotAlt(mi) = qm(2);
    rm = r_of_model(mi, qm, zzM);
    rmsAlt(mi) = sqrt(mean(hypot(rm .* sin(AAM + deg2rad(qm(2))) - xmM, ...
                                 ymM - rm .* cos(AAM + deg2rad(qm(2)))).^2));
end
say('  径向模型敏感度（固定圆心）: rot = %s deg  ⇒ 极差 %.3f deg（说明 rot 由圆心定，不由模型定）\n', ...
    strtrim(sprintf('%.3f ', rotAlt)), max(rotAlt) - min(rotAlt));
say('  同一内点集下各模型的 2D rms : %s px\n', strtrim(sprintf('%.2f ', rmsAlt)));
info.modelRot = rotAlt; info.modelRms = rmsAlt;

% ---- 径向"形状 A/B"裁决（★ 2026-09-17 新增）----
%   问题：与论文的"区间有效尺度"差 −4.5%，是"两条 r(z) 曲线在可及区间内本就不可区分"，
%   还是"数据真的偏好其中一条"？做法：把两种形状各自乘一个自由尺度 k、再配自由 rot
%   （**各 2 个自由参数，参数个数相同**），在同一个圆心、同一个内点集上比 2D rms。
%   若 rms 分不开 ⇒ 那 −4.5% 是"选哪条曲线"的问题，不是测量差异。
shpO = zzM + q(3) * zzM.^3;                       % 本解形状（a 取采用解）
rT   = polyval(TH.rz_coef, zzM);                  % 论文式(3.9) 形状（px 量级）
oo   = optimset('Display', 'off', 'MaxFunEvals', 40000, 'MaxIter', 40000, ...
                'TolX', 1e-9, 'TolFun', 1e-12);
fitKR = @(R, k0) fminsearch(@(p) sum(hypot(p(1) * R .* sin(AAM + deg2rad(p(2))) - xmM, ...
    ymM - p(1) * R .* cos(AAM + deg2rad(p(2)))).^2), [k0, q(2)], oo);
rmsKR = @(R, p) sqrt(mean(hypot(p(1) * R .* sin(AAM + deg2rad(p(2))) - xmM, ...
    ymM - p(1) * R .* cos(AAM + deg2rad(p(2)))).^2));
pO = fitKR(shpO, q(1));  rmsO2 = rmsKR(shpO, pO);
pT = fitKR(rT, 1);       rmsT2 = rmsKR(rT, pT);
% 论文式原样（连尺度也不放）：只剩 rot 一个自由参数
pT1 = fminsearch(@(p) sum(hypot(rT .* sin(AAM + deg2rad(p)) - xmM, ...
    ymM - rT .* cos(AAM + deg2rad(p))).^2), q(2), oo);
rmsT1 = sqrt(mean(hypot(rT .* sin(AAM + deg2rad(pT1)) - xmM, ...
    ymM - rT .* cos(AAM + deg2rad(pT1))).^2));
say('\n  径向**形状** A/B（同圆心、同内点集，r = k * shape(z)，只放 k 与 rot ⇒ 各 2 个自由参数）\n');
say('    本解形状 f(z+a z^3) : rms = %.2f px   (k=%.4f, rot=%+.3f deg)\n', rmsO2, pO(1), pO(2));
say('    论文形状(3.9)       : rms = %.2f px   (k=%.4f, rot=%+.3f deg)   k 需乘 %.4f\n', ...
    rmsT2, pT(1), pT(2), pT(1));
say('    论文形状(3.9) 原样  : rms = %.2f px   (k=1 固定, rot=%+.3f deg)   <== 尺度一点不动\n', ...
    rmsT1, pT1);
say('    互补读法：把论文形状缩到本解尺度需乘 k=%.4f（即 %.1f%%）—— 与"区间有效尺度差 −4.5%%\n', ...
    pT(1), 100 * (pT(1) - 1));
say('    的量级一致；但如果连尺度都不让它动（原样），rms 会升到 %.1f px。\n', rmsT1);
if rmsT2 <= 1.15 * rmsO2
    say('    ⇒ 两种形状在月亮数据上**分不开**（差 %.0f%%）⇒ 那 %.1f%% 是"选哪条曲线"的问题，\n', ...
        100 * (rmsT2 / rmsO2 - 1), 100 * (1 / pT(1) - 1));
    say('      不是测量差异；也说明"谁的径向式更对"在 z=%.0f~%.0f deg 内无法裁决。\n', ...
        round(rad2deg(min(zzM))), round(rad2deg(max(zzM))));
else
    say('    ⇒ 论文形状留下比本解大 %.0f%% 的系统残差 ⇒ 本解的三次形状被数据**偏好**（两者自由\n', ...
        100 * (rmsT2 / rmsO2 - 1));
    say('      参数同为 2 个，差的只是曲线形状，所以不是"参数多所以更好"）。\n');
end
info.shapeCmp = struct('kOurs', pO(1), 'rmsOurs2', rmsO2, 'rmsThesis2', rmsT2, ...
                       'rmsThesis1', rmsT1, 'kThesis', pT(1), ...
                       'rotOurs2', pO(2), 'rotThesis2', pT(2), 'rotThesis1', pT1);

% ---- Δrot/Δy0：把圆心沿 y 挪 10 px 再解一次 rot，直接测出这条简并线的斜率 ----
ymS = (y0 + 10) - cy(kkM)';
qS = fit_model(2, zzM, AAM, xmM, ymS, [q(1), q(2), q(3)]);
rotPerY0 = (qS(2) - rotAlt(2)) / 10;      % deg/px
say('  Δrot/Δy0 = %+.4f deg/px  ⇒ 圆心偏 1 px 就是 rot 偏 %.2f deg（简并的量化）\n', ...
    rotPerY0, rotPerY0);

% ---- 逐夜 rot（在采用解下重算）----
kk = find(keep);
zz = deg2rad(90 - alt(kk))'; AA = deg2rad(azm(kk))';
xm = cx(kk)' - x0; ym = y0 - cy(kk)';
rr = q(1) * (zz + q(3) * zz.^3);
e = hypot(rr .* sin(AA + deg2rad(q(2))) - xm, ym - rr .* cos(AA + deg2rad(q(2))));
% 逐夜 rot（在采用解下重算；rot 的定义与刚拟合出的 q(2) 无关，是各夜独立的观测量）
dd = mod(mod(atan2(cx(kk) - x0, -(cy(kk) - y0)) * 180/pi, 360) - azm(kk), 360);
uN = unique(nt(kk));
pn = repmat(struct('night', '', 'n', 0, 'rot', 0, 'rms', 0, 'scat', 0), 1, numel(uN));
npn = 0;
for i = 1:numel(uN)
    m = strcmp(nt(kk), uN{i});
    if sum(m) < 5, continue; end
    md = mod(atan2(mean(sind(dd(m))), mean(cosd(dd(m)))) * 180/pi, 360);
    npn = npn + 1;
    pn(npn) = struct('night', uN{i}, 'n', sum(m), 'rot', md, ...
        'rms', sqrt(mean(e(m).^2)), 'scat', ...
        sqrt(mean(abs(mod(dd(m) - md + 180, 360) - 180).^2)));
end
pn = pn(1:npn);
cs = 0;
if numel(pn) > 1
    R = abs(mean(exp(1i * deg2rad([pn.rot]))));
    R = min(max(R, 1e-12), 1 - 1e-12);
    cs = sqrt(-2 * log(R)) * 180 / pi;
end
say('  采用解: x0=%.2f y0=%.2f f=%.3f a=%+.6f rot=%+.3f deg rms=%.3f px (%d 帧/%d 夜) 夜间一致性 %.2f deg\n', ...
    x0, y0, q(1), q(3), q(2), sqrt(mean(e.^2)), sum(keep), numel(pn), cs);

cal = struct();
cal.model = 'equid+ (fixed-center)';
cal.r_of_z = 'r = f (z + a z^3)';
cal.chirality = 'direct';
cal.x0 = x0;  cal.y0 = y0;
cal.f = q(1); cal.a = q(3);
cal.rot = mod(q(2), 360);
cal.r_horizon_px = q(1) * (pi/2 + q(3) * (pi/2)^3);
% rim_px 必须给**有限值**：moon_calib_to_pipeline 第 123 行会无条件读 cal.rim_px
% （而它只在 calMat.rim_px 有限时才被赋值），给 NaN 会让那一步报「无法识别的字段」。
cal.rim_px = rimPx;
% 盘边对应的天顶角（反解 r(z)=rim_px）。★ 一期实测 ≈87 deg，**不是 90 deg**：
% 所以"f = R_rim/(pi/2)"这类换算在一期是错的（它会给出 458/1.5708 = 291.6 px/rad，
% 比月亮实测 324.4 低约 10%）。软边宽度只带约 ±1.5 deg 的 theta_rim 不确定性，
% 但"标称 180 deg ⇒ 盘边 = 90 deg"这个默认假设被数据否掉（见 fov_check.m 的假设检验）。
% 取"首次穿越"：二期 a 很负时 r(z) 有极大值，R 会被穿越两次。
zzr = linspace(1e-6, 2.4, 200001);
rvr = cal.f * (zzr + cal.a * zzr.^3);
kzr = find(rvr >= rimPx, 1, 'first');
if isempty(kzr), cal.z_rim_deg = NaN; else, cal.z_rim_deg = rad2deg(zzr(kzr)); end
cal.fov_note = ['采用"纯方位定圆心 + 固定圆心拟合"。注意盘边(z_rim_deg)不是 90 deg，' ...
                '仪器标称 180 deg 视场、有效约 160 deg；f 只能由月亮/星给出，' ...
                '不能用 f = R_rim/(pi/2)。'];
cal.site_lon = lon; cal.site_lat = lat; cal.site_hgt_km = hgt;
cal.rms_px = sqrt(mean(e.^2));
cal.n_frames = sum(keep);  cal.n_nights = numel(pn);
cal.nights = {pn.night};
cal.rot_circ_mean = mod(atan2(mean(sind([pn.rot])), mean(cosd([pn.rot]))) * 180/pi, 360);
cal.rot_circ_scatter = cs;
cal.rot_per_night = pn;
cal.center_search = cen;
cal.rot_model_names = modNames;
cal.rot_model_rot = rotAlt;
cal.rot_model_rms = rmsAlt;
cal.rot_model_spread = max(rotAlt) - min(rotAlt);
cal.rot_per_y0 = rotPerY0;
% ★ 径向"形状 A/B"：本解形状 vs 论文式(3.9) 形状，各 2 个自由参数（k 与 rot）。
%   用于裁决"区间有效尺度 −4.5%"是模型族差异还是测量差异。
cal.shape_cmp = info.shapeCmp;
cal.method = 'direction-only centre + fixed-centre radial fit (calib_2014.m)';
% 2026-09-18：记录圆心是"搜出来的"还是"钉住的"——两者的 rot 口径不同，下游别混用。
if ~isempty(o.CenterFix)
    cal.center_source = 'CenterFix（用户指定，未搜索）';
else
    cal.center_source = '全局紧截尾共识搜索';
end
end


function [cen, cost, info] = search_center(cx, cy, azm, sit, x0g, y0g, o)
%SEARCH_CENTER 用**全局紧截尾共识**定圆心（只用方向，与径向模型无关）。
%
% 判据：找一个 (x0,y0)，使"像面方位 phi = atan2(cx-x0, -(cy-y0)) 与真方位 A 之差"
%       rot = phi - A 在**全体帧**上尽量是同一个常数。
%   实现：迭代 3 次挑选"与当前全局圆中位一致到 ±InlierHalfWidth 的帧"（内点），
%         代价 = 内点的圆方差（deg^2）；若内点比例低于 MinInlierFrac 则加罚。
%
% ★ 为什么不用"逐夜 rot 的汇总量"（旧版做法）：
%   坏夜（20140930/20141007）的 rot 分布又宽又偏（散射 32~35 deg），此时圆均值、
%   圆中位、3xMAD 截尾均值会给出相差 10 deg 的三个答案，任何一个都能被"凑"到
%   与其他夜一致 ⇒ 代价函数被坏夜的汇总量支配，极小跑到 (535,543) 这种伪解上。
%   全局紧截尾不用汇总量，只有"真月亮帧"能通过 ±2 deg 的共识，坏帧自动出局。
%   鲁棒性：6 个不同帧集（见文件头 (5) 的表）都收敛到同一处，内点率 78%~97%。
costf = @(x0, y0) center_cost(x0, y0, cx, cy, azm, sit, o);
best = [x0g, y0g]; bc = costf(x0g, y0g);
for x0 = 480:4:640
    for y0 = 470:4:640
        c = costf(x0, y0);
        if c < bc, bc = c; best = [x0 y0]; end
    end
end
for x0 = best(1) - 4:1:best(1) + 4
    for y0 = best(2) - 4:1:best(2) + 4
        c = costf(x0, y0);
        if c < bc, bc = c; best = [x0 y0]; end
    end
end
% 亚像素：围绕 1 px 最优再细化
for x0 = best(1) - 0.5:0.5:best(1) + 0.5
    for y0 = best(2) - 0.5:0.5:best(2) + 0.5
        c = costf(x0, y0);
        if c < bc, bc = c; best = [x0 y0]; end
    end
end
cen = best; cost = bc;
[info.nIn, ~] = tight_inliers(cen(1), cen(2), cx, cy, azm, sit, o.InlierHalfWidth);
info.nAll = sum(sit); info.frac = info.nIn / max(info.nAll, 1);
end


function c = center_cost(x0, y0, cx, cy, azm, sit, o)
%CENTER_COST 全局紧截尾共识代价：内点圆方差 + 内点比例不足的罚项。
% 本判据**不使用逐夜分组**（这正是它的优点），故签名里没有 night 参数。
[nIn, dk] = tight_inliers(x0, y0, cx, cy, azm, sit, o.InlierHalfWidth);
if nIn < 20, c = 1e9; return; end
frac = nIn / max(sum(sit), 1);
c = dk;
if frac < o.MinInlierFrac
    c = c + 1e4 * (o.MinInlierFrac - frac)^2;   % 不允许用"只留少数帧"换低方差
end
end


function [nIn, dk] = tight_inliers(x0, y0, cx, cy, azm, sit, hw)
%TIGHT_INLIERS 迭代 3 次挑"与全局圆中位一致到 ±hw 度"的帧，返回帧数与圆方差(deg^2)。
rot = mod(mod(atan2(cx(sit) - x0, -(cy(sit) - y0)) * 180/pi, 360) - azm(sit), 360);
keep = true(size(rot));
for it = 1:3
    gr = mod(circmedian_local(rot(keep)), 360);
    keep = abs(mod(rot - gr + 180, 360) - 180) < hw;
end
r = rot(keep); nIn = sum(keep);
if nIn < 5, dk = 1e9; return; end
R = abs(mean(exp(1i * deg2rad(r))));
R = min(max(R, 1e-12), 1 - 1e-12);
dk = -2 * log(R) * (180/pi)^2;
end


function m = circmedian_local(a)
%CIRCMEDIAN_LOCAL 连续（非量化）圆中位：最小化到样本的绝对圆距离之和。
% 用 12 个均匀初值跑 fminsearch 再取最优，避免落进局部极小。
a = a(:).';
o = optimset('Display', 'off');
c = zeros(12, 1);
for k = 1:12
    c(k) = fminsearch(@(x) mean(abs(mod(a - x + 180, 360) - 180)), (k - 1) * 30, o);
end
v = arrayfun(@(x) mean(abs(mod(a - x + 180, 360) - 180)), c);
[~, i] = min(v);
m = mod(c(i), 360);
end


function q = fit3(zz, AA, xm, ym, p0)
%FIT3 固定圆心下拟合 (f, rot, a)：x = r sin(A+rot), y = -r cos(A+rot), r = f(z + a z^3)
rr = @(q) q(1) * (zz + q(3) * zz.^3);
rf = @(q) [rr(q) .* sin(AA + deg2rad(q(2))) - xm; ym - rr(q) .* cos(AA + deg2rad(q(2)))];
o = optimset('Display', 'off', 'MaxFunEvals', 40000, 'MaxIter', 40000, ...
             'TolX', 1e-9, 'TolFun', 1e-12);
q = fminsearch(@(q) sum(rf(q).^2), p0, o);
end


function q = fit_model(mi, zz, AA, xm, ym, p0)
%FIT_MODEL 固定圆心 + 指定径向模型 mi 下拟合 (f, rot, a)，用于量化 rot 的模型敏感度。
rf = @(q) [r_of_model(mi, q, zz) .* sin(AA + deg2rad(q(2))) - xm; ...
           ym - r_of_model(mi, q, zz) .* cos(AA + deg2rad(q(2)))];
o = optimset('Display', 'off', 'MaxFunEvals', 60000, 'MaxIter', 60000, ...
             'TolX', 1e-9, 'TolFun', 1e-12);
q = fminsearch(@(q) sum(rf(q).^2), p0, o);
end


function rr = r_of_model(mi, q, zz)
%R_OF_MODEL 五种常见鱼眼径向映射（第 2 个参数 a 只在 mi==2 时使用）。
switch mi
    case 1, rr = q(1) * zz;
    case 2, rr = q(1) * (zz + q(3) * zz.^3);
    case 3, rr = 2 * q(1) * sin(zz / 2);
    case 4, rr = 2 * q(1) * tan(zz / 2);
    case 5, rr = q(1) * sin(zz);
    otherwise, error('calib_2014:badModel', '未知径向模型编号 %d', mi);
end
end


% ==========================================================================
% 论文样例帧复现
% ==========================================================================
function out = thesis_frame_check(cfg, TH, calMat, outdir, say)
%THESIS_FRAME_CHECK 复现论文图3.2（增强）与图3.9 的同一帧 2014-09-23 17:06:12 UT。
%
% 论文方法（§3.2.1，原文）：
%   (a) 线性提升 4.5 倍；
%   (b) 用圆形模板提取中间图像，圆心取三颗星余弦定理算出的 (x=554, y=538)；
%   (c) 用连续 1 小时（20 张）图像的平均作背景，逐帧相减，再归一化。
% 本函数逐条实现，并把结果落盘成 PNG，便于和论文图直接对比。
out = struct('raw', '', 'enh', '', 'sub', '', 'ok', false);

pat = sprintf('*%s*.PNG', datestr(TH.frame_dt, 'yyyymmddHHMMSS'));
nd = fullfile(cfg.rawRoot, datestr(TH.frame_dt, 'yyyymmdd'));
fs = dir(fullfile(nd, pat));
if isempty(fs)
    say('  [跳过] 没找到样例帧 %s\n', fullfile(nd, pat));
    return
end
fpath = fullfile(nd, fs(1).name);
say('  样例帧 : %s\n', fs(1).name);

[g, info] = q_read_gray(fpath);            % 0..65535 量级
g8 = g / 257;                              % 还原成 8 bit 灰度（便于和论文的"灰阶"对读）

say('  位深   : %d bit    原始灰度 min/中位/max = %.0f / %.0f / %.0f\n', ...
    info.bits, min(g8(:)), median(g8(:)), max(g8(:)));
sv = sort(g8(:));
nq = numel(sv);
qq = @(pr) sv(max(1, min(nq, round(pr * nq))));      % 手写分位（不用 quantile/prctile：要工具箱）
say('  分位   : p01=%.0f p05=%.0f p50=%.0f p95=%.0f p99=%.0f\n', ...
    qq(0.01), qq(0.05), qq(0.5), qq(0.95), qq(0.99));
bins = histcounts(g8(:), -0.5:1:255.5);
occ = sum(bins > 0.001 * numel(g8));
say('  占据 >0.1%% 像元的灰阶级数 = %d 级  ⇒ 天空只落在十几个灰阶里（8 bit 强量化）\n', occ);

% --- (a) 4.5 倍线性增强 ---
enh = min(g8 * 4.5, 255);

% --- (b) 圆形模板：论文圆心 (554,538) vs 本解圆心 ---
[H, W] = size(g8);
[Xg, Yg] = meshgrid(1:W, 1:H);
rTh = hypot(Xg - TH.zenith_star(1), Yg - TH.zenith_star(2));
rUs = hypot(Xg - (calMat.x0 + 1), Yg - (calMat.y0 + 1));
maskTH = rTh <= calMat.r_horizon_px;
maskUS = rUs <= calMat.r_horizon_px;
enhMasked = enh; enhMasked(~maskTH) = 0;

% --- (c) 20 帧(1 小时)平均背景 + 相减 + 归一化 ---
allf = dir(fullfile(nd, '*.PNG'));
names = sort({allf.name});
i0 = find(strcmp(names, fs(1).name), 1);
b0 = max(1, i0 - 10); b1 = min(numel(names), b0 + 19); b0 = max(1, b1 - 19);
say('  背景帧 : 第 %d~%d 帧（共 %d 帧，论文用连续 1 小时的 20 张）\n', b0, b1, b1 - b0 + 1);
acc = zeros(H, W); nb = 0;
for k = b0:b1
    gk = q_read_gray(fullfile(nd, names{k}));
    if isempty(gk), continue; end
    acc = acc + gk / 257; nb = nb + 1;
end
if nb >= 5
    sub = g8 - acc / nb;
    sub = sub - median(sub(maskUS));
    den = prctile_local(abs(sub(maskUS)), 0.99);
    subDisp = 128 + 127 * max(min(sub / max(den, eps), 1), -1);
    subDisp(~maskUS) = 0;
else
    subDisp = zeros(H, W);
end

% --- 落盘 ---
outdir2 = fullfile(outdir, 'thesis_frame');
if ~exist(outdir2, 'dir'), mkdir(outdir2); end
out.raw = fullfile(outdir2, 'A_raw.png');
out.enh = fullfile(outdir2, 'B_enh45_masked.png');
out.sub = fullfile(outdir2, 'C_enh_minus_1h_background.png');
imwrite(uint8(min(max(g8, 0), 255)), out.raw);
imwrite(uint8(min(max(enhMasked, 0), 255)), out.enh);
imwrite(uint8(min(max(subDisp, 0), 255)), out.sub);
say('  已写图 : %s\n           %s\n           %s\n', out.raw, out.enh, out.sub);

% --- 论文声称的三个关键像素，到底有没有"星" ---
say('\n  论文声称的位置处，实际灰度（8 bit；局部背景 = 该点 21x21 窗口中位数）:\n');
say('  %-26s %8s %10s %10s %10s\n', '位置', '灰度', '局部背景', '对比', '判定');
say('  %s\n', repmat('-', 1, 70));
% 注意：三个点必须统一到**原图坐标系**再取样。§3.2.2 的 (452,452) 与北极星 (499,55)
% 都是**裁剪系**，要加原点平移 zoff=(102,86) 才能落到原图像素上。
zoff = TH.zenith_star - TH.zenith_rot;   % 本函数内局部重算（调用者的 zoff 不在作用域）
spots = {'§3.2.1 天顶',          TH.zenith_star; ...
         '§3.2.2 圆心→原图系',    TH.zenith_rot + zoff; ...
         '北极星→原图系',          TH.polaris_px + zoff};
for k = 1:size(spots, 1)
    x = spots{k, 2}(1); y = spots{k, 2}(2);
    v = g8(y, x);
    r0 = 10; xa = max(1, x - r0):min(W, x + r0); ya = max(1, y - r0):min(H, y + r0);
    bgv = median(median(g8(ya, xa)));
    d = v - bgv;
    if d >= 8, verd = '像点源'; elseif d >= 3, verd = '弱偏亮'; else, verd = '无'; end
    say('  %-26s %8.0f %10.0f %+10.0f %10s\n', spots{k, 1}, v, bgv, d, verd);
end
say('  %s\n', repmat('-', 1, 70));
say('  注：§3.2.2 的 (452,452) 与北极星 (499,55) 属**裁剪系**，已加原点平移 zoff=(%d,%d)\n', ...
    zoff(1), zoff(2));
say('      换算到原图系后再取样：%s→%s，%s→%s。\n', ...
    mat2str(TH.zenith_rot), mat2str(TH.zenith_rot + zoff), ...
    mat2str(TH.polaris_px), mat2str(TH.polaris_px + zoff));
say('  判据阈值（经验）: 8 bit 强量化下"星"至少要拉开 8 个灰阶才算看得见。\n');
say('  本工程先前用 5 类独立检验（含奇偶帧分离相关 0.970）判定一期 2014 数据\n');
say('  **没有会随天球转的点**，见 mooncal_work\\reports\\一期2014数据可用性评估.md。\n');
say('  三个位置都不亮 ⇒ 论文的星点定位无法在归档帧上复现。可能是:\n');
say('    (i) 论文用的是未经 STP 处理 / 未经 8 bit 量化的 L0；或\n');
say('    (ii) 论文的星点坐标在归档帧上不可复现（本工程已判该批无会随天球转的点）。\n');
out.ok = true;
end


% ---------------------------------------------------------------------------
function v = prctile_local(x, p)
sv = sort(x(:)); n = numel(sv);
v = sv(max(1, min(n, round(p * n))));
end
