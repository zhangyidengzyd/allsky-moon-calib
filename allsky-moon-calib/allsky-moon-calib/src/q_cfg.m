function cfg = q_cfg(action)
%Q_CFG  ★纯 MATLAB 全链的唯一路径配置★
%
%   cfg = q_cfg();          % 取配置(最常用)
%   q_cfg('show')           % 打印当前配置
%   q_cfg('check')          % 逐个检查路径存在性 + 目录分层体检
%
% 本文件是整条链的"总开关"。换原始图像文件夹、换定标文件、换暗场、换输出盘,
% 都只改这里 —— 五个驱动脚本 (q1~q5) 与 bubble_search / clarity_rank /
% bubble_synth_test 全部从这里取路径, 没有任何别处再写死路径。
%
% 目录约定: rawRoot 下**按日期分子目录**, 帧直接放在日期目录里:
%     <rawRoot>\20250920\ODAZH_DCAI01_AROA_L0_STP_2025092012xxxx_V01.00.PNG
% 各驱动只给**子目录名**(如 '20250920'), 不给全路径。

if nargin < 1 || isempty(action), action = ''; end

% ==========================================================================
% ★★★★★  改路径只改下面这一段  ★★★★★
% ==========================================================================
% 所有路径来自 mc_paths.m（唯一来源）—— 整包不含任何盘符，随目录移动。
% 数据一侧 data/ ，产物一侧 out/ ；详见 README「数据放置」。
MP = mc_paths();
cfg.rawRoot  = MP.raw;        % 原始帧根目录(按日期分子目录)  data/raw/<YYYYMMDD>/*.PNG
cfg.calFile  = MP.calP2Pipe;  % 管线用定标(含 cal 结构: zenith/rot/rz_poly/...)
cfg.calNative= MP.calP2Native;% moon_calib 原生定标(变量 calMat, 0-based)
cfg.darkFile = MP.dark;       % 暗场 PNG。置 '' 会警告, 且耗竭率被 3441 DN 基座稀释一半以上
cfg.projRoot = MP.project;    % 地理投影图输出根目录
cfg.bubRoot  = MP.bubble;     % 泡团处理中间产物(_dep.mat、候选、图)
cfg.logRoot  = MP.logs;       % 文本日志与报告

% 站点与几何约定（同站不改；换站必须改 site，并重跑定标）
cfg.site     = [109.133, 19.526, 0.103];    % [lon_deg, lat_deg, h_km] 官方儋州站(富克 FKT)
cfg.zMax     = 70;                          % 投影限幅天顶角(deg)。两个依据: ①径向标定只在
                                            % z<76 有实测锚点(alt 14~50); ②仪器标称视场180 deg、
                                            % 有效约160 deg ⇒ 可见盘边≈87 deg、可用上限80 deg,
                                            % 70 deg 留足余量。限制来自"分辨率发散+锚点",不来自视场,
                                            % 不要往上抬到75~80。见 fov_check.m / 报告 §5
cfg.resDeg   = 0.02;                        % 画布分辨率(deg) ≈ 2.22 km
cfg.hShell   = 250;                         % 气辉层高度(km)
cfg.pattern  = 'ODAZH_DCAI01_AROA_L0_STP_*.PNG';      % 帧文件通配
cfg.dtSec    = 180;                         % 采样间隔(s)，用于时长/漂移率换算
% ==========================================================================

switch lower(char(action))
    case ''
        return
    case 'show'
        f = fieldnames(cfg);
        for k = 1:numel(f)
            v = cfg.(f{k});
            if ischar(v)
                fprintf('%-11s = %s\n', f{k}, v);
            elseif isscalar(v) && isnumeric(v)
                fprintf('%-11s = %g\n', f{k}, v);
            else
                fprintf('%-11s = %s\n', f{k}, mat2str(v, 6));
            end
        end
    case 'check'
        q_check(cfg);
    otherwise
        error('q_cfg:badAction', '只支持 ''show'' 或 ''check''。');
end
end

% ==========================================================================
function q_check(cfg)
% 逐个路径做存在性检查，并把"目录分层是否对"一并体检。
P = @(varargin) fprintf(varargin{:});
P('\n===== q_cfg 路径体检 =====\n');
chk(cfg.calFile,   'file', '管线定标');
chk(cfg.calNative, 'file', '原生定标');
if isempty(cfg.darkFile)
    P('  [警告] 暗场为空: 耗竭率会被 3441 DN 基座稀释一半以上, 结果只能看形态不能看振幅。\n');
else
    chk(cfg.darkFile, 'file', '暗场');
end
chk(cfg.projRoot, 'new', '投影输出');
chk(cfg.bubRoot,  'new', '泡团输出');
chk(cfg.logRoot,  'new', '日志目录');
chk(cfg.rawRoot,  'dir', '原始帧根目录');

% 分层约定体检：rawRoot 下应当是"日期目录"，日期目录里应当直接有帧
if exist(cfg.rawRoot, 'dir')
    d = dir(cfg.rawRoot);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    if isempty(d)
        P('  [警告] rawRoot 下没有任何子目录 —— 期望"按日期分子目录"。\n');
    else
        n = 0; nested = 0;
        for k = 1:numel(d)
            f = dir(fullfile(cfg.rawRoot, d(k).name, cfg.pattern));
            if isempty(f), f = dir(fullfile(cfg.rawRoot, d(k).name, '*.PNG')); end
            n = n + numel(f);
            s = dir(fullfile(cfg.rawRoot, d(k).name));
            s = s([s.isdir] & ~ismember({s.name}, {'.', '..'}));
            if ~isempty(s) && isempty(f), nested = nested + 1; end
        end
        P('  rawRoot 下 %d 个日期目录, 共 %d 帧。\n', numel(d), n);
        if nested > 0
            P('  [警告] %d 个日期目录里没有直接放帧、而是又套了一层目录。\n', nested);
            P('         本链要求: <rawRoot>\\<日期>\\*.PNG\n');
        end
    end
end

    function chk(p, kind, name)
        switch kind
            case 'file', ok = (exist(p, 'file') == 2);
            case 'dir',  ok = (exist(p, 'dir') == 7);
            case 'new',  ok = (exist(p, 'dir') == 7);
        end
        if ok
            P('  [OK]   %-14s %s\n', name, p);
        elseif strcmp(kind, 'new')
            P('  [待建] %-14s %s   (首次运行时自动建)\n', name, p);
        else
            P('  [缺!]  %-14s %s\n', name, p);
        end
    end
end
