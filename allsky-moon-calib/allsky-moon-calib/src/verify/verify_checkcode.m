function acc = verify_checkcode(verbose)
%VERIFY_CHECKCODE  对包内每个 .m 跑一次 MATLAB 静态检查（checkcode）。
%
%   verify_checkcode          % 逐文件打印
%   verify_checkcode(false)   % 只打印有问题的
%   acc = verify_checkcode()  % acc.nFile / acc.nWarn / acc.bad / acc.legacy
%
% 判据
%   非 LEGACY 文件里出现任何 checkcode 条目 ⇒ 计入 acc.bad（"有问题的文件"数）。
%   moon_calib.m 是 94 KB 的遗留大脚本，checkcode 会给一批"预分配内存 / 未用输入参数"
%   之类**性能与风格提示**（不是错误，也不影响数值），因此单列进 LEGACY，
%   只报条数、不判失败。
%
% 为什么单列一个文件：q4_verify 的 ② 段与本函数是同一件事；抽出来之后
%   `run_demo('verify')` 可以在**没有原始帧**的情况下也跑静态检查。

if nargin < 1 || isempty(verbose), verbose = true; end

LEGACY = {'moon_calib.m'};

[kit, kitF] = package_file_list();

acc.nFile = numel(kit);
acc.nWarn = 0;
acc.bad   = 0;
acc.legacy = struct('name', {}, 'n', {});

fprintf('\n########## 静态检查 checkcode（%d 个文件）##########\n', acc.nFile);
for k = 1:numel(kit)
    f = kitF{k};
    if exist(f, 'file') ~= 2
        fprintf('  [缺] %s\n', kit{k});
        acc.bad = acc.bad + 1;
        continue
    end
    r = checkcode(f, '-struct');
    isLegacy = any(strcmp(kit{k}, LEGACY));
    if isempty(r)
        if verbose, fprintf('  [OK]   %s\n', kit{k}); end
    elseif isLegacy
        acc.legacy(end + 1) = struct('name', kit{k}, 'n', numel(r));
        fprintf('  [信息] %s —— %d 条，均为遗留脚本的性能/风格提示，不参与判读\n', ...
                kit{k}, numel(r));
        if verbose
            for q = 1:numel(r)
                fprintf('         L%d: %s\n', r(q).line, r(q).message);
            end
        end
    else
        acc.nWarn = acc.nWarn + numel(r);
        acc.bad   = acc.bad + 1;
        fprintf('  [%d 条] %s\n', numel(r), kit{k});
        for q = 1:numel(r)
            fprintf('         L%d: %s\n', r(q).line, r(q).message);
        end
    end
end

fprintf('  ---- 合计：%d 个文件，非遗留告警 %d 条，有问题文件 %d 个\n', ...
        acc.nFile, acc.nWarn, acc.bad);
if acc.bad == 0
    fprintf('  [OK]   全部通过（moon_calib.m 按遗留脚本单列）。\n');
end
end
