function [L, F] = package_file_list()
%PACKAGE_FILE_LIST  包内全部 .m 的名称与完整路径（q4_verify 的静态检查用）
%
%   [L, F] = package_file_list()
%     L : N x 1 文件名（不含路径），已排序
%     F : N x 1 完整路径
%
% 为什么单列一个函数：静态检查的清单最容易与实际文件集脱节，一旦脱节就会报
% "[缺] xxx.m"，看起来像链路坏了、其实只是清单过时。这里按本文件所在目录**实时枚举**
% src/ 与 src/verify/ 两处，永远与包同步。

here = fileparts(mfilename('fullpath'));       % <pkg>/src
F = {};
d1 = dir(fullfile(here, '*.m'));
d2 = dir(fullfile(here, 'verify', '*.m'));
for d = {d1, d2}
    dd = d{1};
    for k = 1:numel(dd)
        if strcmp(dd(k).name, 'package_file_list.m'), continue; end
        F{end + 1, 1} = fullfile(dd(k).folder, dd(k).name);   %#ok<AGROW>
    end
end
F = sort(F);
[~, n, ~] = cellfun(@(p) fileparts(p), F, 'UniformOutput', false);
L = n(:);
end
