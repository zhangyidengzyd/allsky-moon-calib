function P = mc_paths()
%MC_PATHS  ★包内所有数据与产物路径的唯一来源★
%
%   P = mc_paths()
%
% 为什么要它
%   整包必须能"拷到任意位置、不改一行路径就能跑"。所有路径都从 pkg_root()
%   派生，因此**任何文件里都不许再出现盘符**。要改目录结构，只改这一个文件。
%
% 约定
%   data/raw/<YYYYMMDD>/*.PNG   原始帧（未随包提供，见 README「数据放置」）
%   data/dark/                  暗场 PNG
%   data/calib/                 定标结果与逐帧月亮量测缓存（随包提供的参考解）
%   out/                        所有可再生的产物（重跑会覆盖）
%
% ★ 重跑 start_calib / calib_2014 会**覆盖** data/calib/ 下的参考解。
%   需要保留请先整个拷一份。

r = pkg_root();
P.root = r;
P.src  = fullfile(r, 'src');
P.data = fullfile(r, 'data');
P.docs = fullfile(r, 'docs');

% ---- 数据一侧（输入）----
P.raw   = fullfile(P.data, 'raw');
P.dark  = fullfile(P.data, 'dark', 'ODAZH_DCAI01_DFOA_AUX_STP_20260516115948_V01.00.PNG');
P.calib = fullfile(P.data, 'calib');
P.calP2Native  = fullfile(P.calib, 'phase2_native.mat');      % 二期：calMat，0-based
P.calP2Pipe    = fullfile(P.calib, 'phase2_pipe.mat');        % 二期：cal 结构，管线用
P.measP2       = fullfile(P.calib, 'phase2_meas.mat');        % 二期：逐帧月核量测缓存
P.calP14Native = fullfile(P.calib, 'phase1_2014_native.mat'); % 一期 2014
P.calP14Pipe   = fullfile(P.calib, 'phase1_2014_pipe.mat');
P.measP14      = fullfile(P.calib, 'phase1_2014_meas.mat');
P.calP13Native = fullfile(P.calib, 'phase1_2013_native.mat'); % 一期 2013
P.calP13Pipe   = fullfile(P.calib, 'phase1_2013_pipe.mat');
P.measP13      = fullfile(P.calib, 'phase1_2013_meas.mat');

% ---- 产物一侧（输出，可再生）----
P.out        = fullfile(r, 'out');
P.project    = fullfile(P.out, 'project');
P.bubble     = fullfile(P.out, 'bubble');
P.logs       = fullfile(P.out, 'logs');
P.mooncal    = fullfile(P.out, 'mooncal');                     % moon_calib 的产物根
P.mooncalOut = fullfile(P.mooncal, 'mooncal_out_matlab');
end
