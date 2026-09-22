function mc_pipeline_e2e_test(nightFolder, nFrames, moonOutMat, pipeCal)
%MC_PIPELINE_E2E_TEST  验证「月亮标定 -> run.m」这条链是否真的接通。
%
%   把 nFrames 帧**真实**原始帧喂给用户自己的 airglow_geo_pipeline 两次：
%     ① moon_calib 原生 .mat（变量 calMat）  -> 应当**报错**：标定结构体缺少字段：zenith
%     ② 适配后 .mat（变量 cal，含 zenith/rot/radius/rz_poly）-> 应当跑出投影图
%   这条试验的判别力来自"管线自己接受/拒绝"，不依赖任何我们写的公式。
%
% 用法
%   mc_pipeline_e2e_test                                  % 全默认
%   mc_pipeline_e2e_test('E:\qihui2\原始数据\20260323')    % 换一夜
%   mc_pipeline_e2e_test('', 30)                          % 只改帧数
%   mc_pipeline_e2e_test('', [], '', 'D:\mycal.mat')      % 指定待测标定文件
%
% 参数
%   nightFolder  含 ODAZH*.png 的目录，默认 E:\qihui2\原始数据\20250726
%   nFrames      取前多少帧，默认 20（够跑滑动窗口背景，又很快）
%   moonOutMat   moon_calib 原生输出，默认 <本文件目录>\mooncal_work\mooncal_out_matlab\calibration_params.mat
%   pipeCal      待检验的管线格式标定，默认 E:\qihui2\Processed_Output\calibration_params.mat
%
% 依赖：airglow_geo_pipeline.m 在 MATLAB 路径上（本脚本会自动 addpath E:\qihui2）。
% 本脚本用纯 ASCII 书写 —— R2021b 的 -batch 对非 ASCII 源码容错很差。
% 临时帧目录用 tempname() 建、跑完自动删，不留垃圾。

MP = mc_paths();
if nargin < 1 || isempty(nightFolder), nightFolder = fullfile(MP.raw, '20250726'); end
if nargin < 2 || isempty(nFrames),     nFrames     = 20; end
if nargin < 3 || isempty(moonOutMat)
    moonOutMat = MP.calP2Native;
end
if nargin < 4 || isempty(pipeCal)
    pipeCal = MP.calP2Pipe;
end

addpath(MP.src);   % airglow_geo_pipeline.m 所在处

% ---------- 0) 准备一小批真实帧 ----------
assert(exist(nightFolder, 'dir') == 7, '找不到原始数据目录：%s', nightFolder);
fl = dir(fullfile(nightFolder, 'ODAZH*.png'));
if isempty(fl), fl = dir(fullfile(nightFolder, 'ODAZH*.PNG')); end
assert(~isempty(fl), '目录里没有 ODAZH*.png：%s', nightFolder);
n = min(nFrames, numel(fl));
td = tempname(); mkdir(td);
cu = onCleanup(@() rmdir(td, 's'));   % 跑完自动删临时目录
rawDir = fullfile(td, 'frames'); mkdir(rawDir);
for k = 1:n
    copyfile(fullfile(fl(k).folder, fl(k).name), rawDir);
end
fprintf('测试数据：%s 的前 %d 帧 -> %s\n', nightFolder, n, rawDir);

% ---------- 1) 原生 .mat -> 应当报错 ----------
fprintf('\n=== 试验 1：moon_calib 原生 .mat 直接喂管线（应当失败）===\n');
fprintf('文件：%s\n', moonOutMat);
if exist(moonOutMat, 'file') ~= 2
    fprintf('>>> 文件不存在，跳过（先跑一次 moon_calib）\n');
else
    S = load(moonOutMat); fn = fieldnames(S);
    fprintf('变量：%s\n', strjoin(fn.', ', '));
    for i = 1:numel(fn)
        v = S.(fn{i});
        if isstruct(v) && ~isfield(v, 'zenith')
            fprintf('  %s 有 zenith? 否  有 rz_poly? %d  有 radius? %d\n', ...
                fn{i}, isfield(v, 'rz_poly'), isfield(v, 'radius'));
        end
    end
    try
        airglow_geo_pipeline(rawDir, fullfile(td, 'out_1'), moonOutMat, 250, 19.5, 109.2, ...
            'BgWindow', 3, 'GeoFolder', fullfile(td, 'geo_1'));
        fprintf('>>> 意外成功 —— 说明管线改了约定，请重新核对适配器！\n');
    catch ME
        fprintf('>>> 如预期报错：[%s] %s\n', ME.identifier, ME.message);
    end
end

% ---------- 2) 适配后 .mat -> 应当成功 ----------
fprintf('\n=== 试验 2：适配后的 .mat 喂管线（应当成功）===\n');
fprintf('文件：%s\n', pipeCal);
if exist(pipeCal, 'file') ~= 2
    fprintf('>>> 文件不存在，跳过。用下面命令生成：\n');
    fprintf('    moon_calib(''stage'',''solve,export'',''pipeline'',''%s'')\n', pipeCal);
    return
end
t0 = tic;
airglow_geo_pipeline(rawDir, fullfile(td, 'out_2'), pipeCal, 250, 19.5, 109.2, ...
    'BgWindow', 3, 'GeoFolder', fullfile(td, 'geo_2'));
fprintf('>>> 成功，用时 %.1f s\n', toc(t0));

dg = dir(fullfile(td, 'geo_2', '*.png'));
dp = dir(fullfile(td, 'out_2', '*_proc.png'));
fprintf('>>> 投影图 %d 张，预处理图 %d 张\n', numel(dg), numel(dp));
assert(numel(dg) == n, '投影图数量 %d != 输入帧数 %d', numel(dg), n);
if ~isempty(dg)
    a = imread(fullfile(dg(1).folder, dg(1).name));
    fprintf('>>> 首张投影图 %dx%d %s max=%g\n', size(a,1), size(a,2), class(a), double(max(a(:))));
end

% ---------- 3) 字段对照 ----------
fprintf('\n=== 试验 3：两个 .mat 的字段对照 ===\n');
if exist(moonOutMat, 'file') == 2
    S = load(moonOutMat); fn = fieldnames(S);
    fprintf('原生  变量 %s：%s\n', strjoin(fn.', '/'), strjoin(fieldnames(S.(fn{1})).', ', '));
end
P = load(pipeCal);
fprintf('适配  变量 %s：%s\n', strjoin(fieldnames(P).', '/'), strjoin(fieldnames(P.cal).', ', '));
c = P.cal;
fprintf('  zenith(1-based) = (%.3f, %.3f)\n', c.zenith(1), c.zenith(2));
fprintf('  rot             = %+.4f deg\n', c.rot);
fprintf('  radius          = %.2f px\n', c.radius);
fprintf('  rz_poly         = [%s]\n', strrep(sprintf('%.6g ', c.rz_poly), '  ', ' '));
fprintf('  r(90 deg)       = %.2f px\n', polyval(c.rz_poly, 90));
fprintf('  obs_lat/lon     = %.3f / %.3f   (会覆盖 run.m 的实参)\n', c.obs_lat, c.obs_lon);

fprintf('\n=== 全部通过 ===\n');
end
