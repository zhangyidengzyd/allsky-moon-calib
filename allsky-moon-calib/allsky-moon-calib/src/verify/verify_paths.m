function acc = verify_paths()
%VERIFY_PATHS  包的三项硬自检：① 可执行代码里没有盘符 ② 依赖零缺口 ③ 数据文件齐备。
%
%   verify_paths
%   acc = verify_paths()      % acc.n / acc.bad
%
% ① 判别口径：**可执行代码**里的 `盘符:\`。注释行里的旧路径**不算**（保留作溯源），
%    string 字面量里的也不算（本文件自己的正则就写在字符串里）。
% ② 判据：包内每个 .m 里出现的 `名字(` 形式的调用，若既不是包内定义的函数、
%    MATLAB 又不认识（exist == 0），就算缺口。函数形参与赋值左侧会被先剔除，
%   所以不会把普通变量误报成函数。
% ③ 判据：mc_paths() 里所有"输入侧"路径存在（产物侧不存在是正常的，会自建）。
%
% ★ 设计上刻意不依赖任何工具箱，纯 base MATLAB，便于审稿人在任何 R2018b+ 上跑。

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % <pkg>
addpath(fullfile(root, 'src'));
MP = mc_paths();

acc.n = 0; acc.bad = 0;

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  包自检  包根 %s\n', root);
fprintf('%s\n', repmat('=', 1, 78));

% 待检文件：src/*.m + src/verify/*.m + 包根的 run_demo.m
F = [dir(fullfile(MP.src, '*.m')); dir(fullfile(MP.src, 'verify', '*.m')); ...
     dir(fullfile(root, '*.m'))];
F = F(~cellfun(@isempty, {F.name}));

% ==========================================================================
% ① 可执行代码里的盘符
% ==========================================================================
fprintf('\n--- ① 绝对路径（盘符）扫描：可执行代码\n');
hits = {};
defs = {};        % 包内定义的全部函数名
for k = 1:numel(F)
    p = fullfile(F(k).folder, F(k).name);
    L = readlines_local(p);
    for i = 1:numel(L)
        code = strip_code(L{i});
        if ~isempty(regexp(code, '[A-Za-z]:\\', 'once'))
            hits{end + 1} = sprintf('%s:%d', F(k).name, i);   %#ok<AGROW>
            fprintf('  [!!] %s:%d  %s\n', F(k).name, i, strtrim(L{i}));
        end
    end
    if ~isempty(F(k).name)
        nf = local_func_names(L);
        defs = [defs(:); nf(:)];   %#ok<AGROW>
    end
end
acc.n = acc.n + 1;
if isempty(hits)
    fprintf('  [OK] 可执行代码含盘符的行：0（注释行不计）\n');
else
    fprintf('  合计 %d 行。\n', numel(hits)); acc.bad = acc.bad + 1;
end

% ==========================================================================
% ② 依赖零缺口
% ==========================================================================
fprintf('\n--- ② 依赖闭包：包内调用是否都能解析\n');
unresolved = {};
for k = 1:numel(F)
    p = fullfile(F(k).folder, F(k).name);
    L = readlines_local(p);
    txt = strjoin(L, sprintf('\n'));
    code = strip_code(txt);
    % 变量集合：函数形参/返回值 + 赋值左侧 + for 循环变量
    vars = local_var_names(L);
    calls = regexp(code, '(?<![\w.])([A-Za-z]\w*)\s*\(', 'tokens');
    calls = unique(cellfun(@(c) c{1}, calls, 'UniformOutput', false));
    for i = 1:numel(calls)
        nm = calls{i};
        if any(strcmp(nm, vars)), continue; end
        if any(strcmp(nm, defs)), continue; end
        if exist(nm) ~= 0, continue; end          % MATLAB 自带（函数/builtin/类）
        unresolved{end + 1} = sprintf('%s  <- %s', nm, F(k).name);   %#ok<AGROW>
    end
end
acc.n = acc.n + 1;
if isempty(unresolved)
    fprintf('  [OK] 未解析的调用：0（包内 %d 个函数 + MATLAB 自带全部可解析）\n', numel(unique(defs)));
else
    fprintf('  [!!] %d 个调用无法解析：\n', numel(unresolved));
    for i = 1:numel(unresolved), fprintf('        %s\n', unresolved{i}); end
    acc.bad = acc.bad + 1;
end

% ==========================================================================
% ③ 数据文件齐备（输入侧）
% ==========================================================================
fprintf('\n--- ③ 输入侧数据是否齐备\n');
need = { ...
  '暗场',                 MP.dark; ...
  '二期 原生定标',        MP.calP2Native; ...
  '二期 管线定标',        MP.calP2Pipe; ...
  '二期 月核量测缓存',    MP.measP2; ...
  '一期2014 原生定标',    MP.calP14Native; ...
  '一期2014 管线定标',    MP.calP14Pipe; ...
  '一期2014 月核量测缓存',MP.measP14; ...
  '一期2013 原生定标',    MP.calP13Native; ...
  '一期2013 管线定标',    MP.calP13Pipe; ...
  '一期2013 月核量测缓存',MP.measP13};
for k = 1:size(need, 1)
    acc.n = acc.n + 1;
    if exist(need{k, 2}, 'file') == 2
        fprintf('  [OK]   %-22s %8.1f KB\n', need{k, 1}, dir(need{k, 2}).bytes / 1024);
    else
        fprintf('  [缺!]  %-22s %s\n', need{k, 1}, need{k, 2}); acc.bad = acc.bad + 1;
    end
end
acc.n = acc.n + 1;
if exist(MP.raw, 'dir') == 7
    d = dir(MP.raw);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    nfr = 0;
    for k = 1:numel(d)
        f = dir(fullfile(MP.raw, d(k).name, '*.PNG'));
        nfr = nfr + numel(f);
    end
    fprintf('  [OK]   原始帧目录 %s ：%d 个日期目录 / %d 帧\n', rel(MP.root, MP.raw), numel(d), nfr);
    if numel(d) == 0
        fprintf('         （空是正常的——原始帧未随包提供，见 README「数据放置」）\n');
    end
else
    fprintf('  [缺!]  原始帧目录 %s\n', MP.raw); acc.bad = acc.bad + 1;
end

% ==========================================================================
fprintf('\n%s\n', repmat('=', 1, 78));
if acc.bad == 0
    fprintf('  三项自检全部通过（%d 项）。\n', acc.n);
else
    fprintf('  %d 项里有 %d 项未通过。\n', acc.n, acc.bad);
end
fprintf('%s\n\n', repmat('=', 1, 78));
end

% ==========================================================================
%  以下为局部工具
% ==========================================================================
function L = readlines_local(p)
fid = fopen(p, 'r', 'n', 'UTF-8');
if fid < 0, L = {}; return; end
c = onCleanup(@() fclose(fid));
L = {};
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    L{end + 1} = ln;   %#ok<AGROW>
end
end

function out = strip_code(s)
%STRIP_CODE 去掉注释与字符串字面量，只留可执行代码。
% ' 的歧义（转置 vs 字符串）用启发式：前一非空字符是字母/数字/)/]/}/. 时按转置处理。
out = ''; i = 1; n = numel(s); inStr = false;
prevSig = '';      % 前一个非空字符
while i <= n
    c = s(i);
    if inStr
        if c == ''''
            if i + 1 <= n && s(i + 1) == ''''
                i = i + 2; continue
            end
            inStr = false;
        end
        i = i + 1; continue
    end
    if c == '%'
        break                                   % 注释：本行余下全丢
    end
    if c == ''''
        if ~isempty(prevSig) && ...
           (isletter(prevSig) || (prevSig >= '0' && prevSig <= '9') || ...
            any(prevSig == ')]}.'))
            out(end + 1) = c;                   %#ok<AGROW>  % 转置符，保留
            i = i + 1; continue
        end
        inStr = true; i = i + 1; continue
    end
    out(end + 1) = c;                           %#ok<AGROW>
    if ~isspace(c), prevSig = c; end
    i = i + 1;
end
end

function names = local_func_names(L)
%LOCAL_FUNC_NAMES 从各行里抽出 `function` 定义的名字（顶层 + 局部）。
names = {};
for i = 1:numel(L)
    s = strtrim(L{i});
    if ~strncmp(s, 'function', 8), continue; end
    t = regexp(s, '^\s*function\s+(?:\[[^\]]*\]|[\w~]+(?:\s*,\s*[\w~]+)*)?\s*=?\s*(\w+)\s*\(', 'tokens', 'once');
    if ~isempty(t), names{end + 1} = t{1}; end   %#ok<AGROW>
end
end

function vars = local_var_names(L)
%LOCAL_VAR_NAMES 收集"明显是变量"的名字：函数形参/返回值、赋值左侧、for 变量。
vars = {};
for i = 1:numel(L)
    s = strtrim(L{i});
    if strncmp(s, '%', 1), continue; end
    s = strip_code(s);
    if isempty(strtrim(s)), continue; end
    % function [a,b] = f(x,y)
    t = regexp(s, '^\s*function\s+(\[[^\]]*\]|\w+)?\s*=?\s*(\w+)\s*\(([^)]*)\)', 'tokens', 'once');
    if ~isempty(t)
        for part = {t{1}, t{3}}
            if ~isempty(part{1})
                w = regexp(part{1}, '[A-Za-z]\w*', 'match');
                vars = [vars, w];   %#ok<AGROW>
            end
        end
        continue
    end
    % for i = ... / parfor
    t = regexp(s, '^\s*(?:par)?for\s+(\w+)\s*=', 'tokens', 'once');
    if ~isempty(t), vars{end + 1} = t{1}; continue; end          %#ok<AGROW>
    % 赋值左侧：[a, b] = ...  或  a = ...  （不是 ==，也不是 >=/<=/~=）
    t = regexp(s, '^\s*(\[[^\]]*\]|[A-Za-z]\w*(?:\{\d+\})?(?:\([^)]*\))?)\s*=(?!=)', 'tokens', 'once');
    if ~isempty(t)
        vars = [vars, regexp(t{1}, '[A-Za-z]\w*', 'match')];     %#ok<AGROW>
    end
end
vars = unique(vars);
end

function p = rel(root, f)
p = strrep(f, [root filesep], '');
end
