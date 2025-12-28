clear;clc;
% 定义参数
k = 4; % 假设 k=1，您可以根据需要调整

% 定义范围
t = linspace(-0.5, 0.5, 1000); % 时间，200个点
f = linspace(0, 100, 1000); % 频率 f，200个点
f_optloc = 70; % 频率 f_optloc，200个点

% 方法1：先创建空数组再赋值（最清晰）
h = zeros( length(f), length(t)); % 注意维度顺序

% 三层循环填充数据（确保每个点都计算）
    for j = 1:length(f)
        for k_idx = 1:length(t)
            % 计算 delta_f
            delta_f = f_optloc - f(j);
            
            % 计算 sigma
            sigma = (f_optloc + abs(delta_f))/k;%sigma=sigma^-1
            
            % 计算 h 值（确保 sigma > 0 避免除零）
            if sigma > 0
                h( j, k_idx) = 1/sqrt(2*pi).*sigma * exp(-t(k_idx)^2/2.*sigma^2);
            else
                h( j, k_idx) = 0; % 或 NaN
            end
        end
    end

% 提取切片数据
%slice_index_f_optloc = find(abs(f_optloc-80) == min(abs(f_optloc-80)));
%slice_f = squeeze( :, :));

% 创建网格
[T, F] = meshgrid(t, f);

% 绘制三维曲面图
figure('Position', [100, 100, 800, 600]);
surf(T, F, h, 'EdgeColor', 'none');
xlabel('时间 t (s)', 'FontSize', 12);
ylabel('f_{optloc} (Hz)', 'FontSize', 12);
zlabel('h(t, f, f_{optloc})', 'FontSize', 12);
title('三维曲面: h(t, f-inst= 80 Hz, f)');

% 美化图形
colormap('jet');
colorbar;
shading interp; % 平滑着色
lighting gouraud;
light('Position', [1, 1, 1]);
light('Position', [-1, -1, -1]);
material dull; % 或 shiny, metal
grid on;
view(45, 30); % 设置视角
axis tight;
