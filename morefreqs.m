% 设置音频文件所在的文件夹路径
myFolder = 'E:\dsp大作业\data\空弦音'; % 替换为您音频文件夹的实际路径

% 检查文件夹是否确实存在
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% 获取文件夹中所有.wav文件的列表
filePattern = fullfile(myFolder, '*.wav');
wavFiles = dir(filePattern);

%设置音频文件保存路径
savepath = 'E:\dsp大作业\picture\空弦音';
savepath1 = 'E:\dsp大作业\picture\空弦音\时域波形和频谱图';
%%
% %循环读取文件夹中的所有.wav文件，绘制时域波形图和频谱图
for k = 1:length(wavFiles)
    baseFileName = wavFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    %读取音频文件
    [y,fs] = audioread(fullFileName);
    %绘制时域波形图
    figure(1);
    subplot(2,1,1);
    plot(y);
    title('时域波形图');
    xlabel('时间');
    ylabel('幅度');
    %绘制频谱图
    subplot(2,1,2);
    %做傅里叶变换，注意这里的y是复数，所以要取绝对值并且只取前一半
    Y = abs(fft(y));
    Y = Y(1:length(Y)/2);
    plot(Y);
    title('频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0 10000]);
    %保存图片
    saveas(gcf,fullfile(savepath1,baseFileName(1:end-4)),'png');
end

%%
% %读取音频文件r1a.wav
% filename = 'E:\dsp大作业\data\空弦音\r1a.wav';

% %读取音频文件
% [y,fs] = audioread(filename);

% %做傅里叶变换，注意这里的y是复数，所以要取绝对值并且只取前一半
% Y = abs(fft(y));
% Y = Y(1:length(Y)/2);

% % %绘制频谱图
% % figure(1);
% % plot(Y);
% % title('频谱图');
% % xlabel('频率');
% % ylabel('幅度');

% %freq为频频率，siganl为对应的频域幅度.所有频率的幅度都在signal中
% siganl = Y;
% freq = (0:length(Y)-1)*fs/length(Y);
% %统计signal中的非零元素个数
% count1 = 0;
% for i = 1:length(siganl)
%     if siganl(i) ~= 0
%         count1 = count1 + 1;
%     end
% end


% %绘制频谱图
% figure(2);
% plot(freq,siganl);
% title('频谱图');
% xlabel('频率');
% ylabel('幅度');
% xlim([0 5000]);

% %多倍频法定位凸峰
% new_signal = processSignal(siganl, freq);


% %绘制new_signal的频谱图
% figure(3);
% plot(freq,new_signal);
% title('频谱图');
% xlabel('频率');
% ylabel('幅度');
% %调整x范围为0-1000
% xlim([0 5000]);
% %new_freq表示new_signal中非零元素对应的频率
% %new_H表示new_signal中非零元素对应的幅度
% new_freq = [];
% new_H = [];
% %统计new_signal中的非零元素个数
% count2 = 0;
% for i = 1:length(new_signal)
%     if new_signal(i) ~= 0
%         count2 = count2 + 1;
%         new_freq(count2) = freq(i);
%         new_H(count2) = new_signal(i);
%     end
% end
% % %打印非零元素个数
% % printf("非零元素个数为：%d\n",count);
%%

%循环读取文件夹中的所有.wav文件，仿照上面的代码计算频谱图，并画出处理前后的频谱图，且把全部处理后的频率特性保存到一个数组中
%这里的数组是一个二维数组，第一维表示文件序号，第二维表示频率序号，数组元素表示对应的幅度
new_freq = [];
new_H = [];
for k = 1:length(wavFiles)
    baseFileName = wavFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    %读取音频文件
    [y,fs] = audioread(fullFileName);

    %做傅里叶变换，注意这里的y是复数，所以要取绝对值并且只取前一半
    Y = abs(fft(y));
    Y = Y(1:length(Y)/2);

    % %绘制频谱图
    % figure(1);
    % plot(Y);
    % title('频谱图');
    % xlabel('频率');
    % ylabel('幅度');

    %freq为频频率，siganl为对应的频域幅度.所有频率的幅度都在signal中
    siganl = Y;
    freq = (0:length(Y)-1)*fs/length(Y);

    %多倍频法定位凸峰
    new_signal = processSignal(siganl, freq);

    %绘制siganl的频谱图和new_signal的频谱图
    figure(2);
    subplot(2,1,1);
    plot(freq,siganl);
    title('频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0 5000]);
    subplot(2,1,2);
    plot(freq,new_signal);
    title('频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0 5000]);
    %保存图片
    saveas(gcf,fullfile(savepath,baseFileName(1:end-4)),'png');

    %把数据保存到数组中
    %new_freq[k]保存的是第k个文件的频率序列
    %new_H[k]保存的是第k个文件的幅度序列
    new_freq{k} = freq;
    new_H{k} = new_signal;
     
end
%%
%把new_freq和new_H保存当前目录的kong.xlxs文件中
%保存
% %使用writematrix函数
% writematrix(new_freq,'kong.xlsx','Sheet',1);
% writematrix(new_H,'kong.xlsx','Sheet',2);
%改用writecell函数
%只写入new_H中的非零元素
% writecell(new_freq,'kong.xlsx','Sheet',1);
% writecell(new_H,'kong.xlsx','Sheet',2);
%%
%把new_H的1,2,3取平均值，然后画出频谱图名字为r1
%把new_H的4,5,6取平均值，然后画出频谱图名字为r4
%把new_H的7,8,9取平均值，然后画出频谱图名字为r6
%补0把维度调成一样的
max_length = 0;
%找到最长的一列
for i = 1:length(new_H)
    if length(new_H{i}) > max_length
        max_length = length(new_H{i});
    end
end
%把所有的列都补0
for i = 1:length(new_H)
    for j = length(new_H{i})+1:max_length
        new_H{i}(j) = 0;
    end
end
%把new_freq也补0
for i = 1:length(new_freq)
    for j = length(new_freq{i})+1:max_length
        new_freq{i}(j) = 0;
    end
end

r1_H = (new_H{1} + new_H{2} + new_H{3})/3;
r4_H = (new_H{4} + new_H{5} + new_H{6})/3;
r6_H = (new_H{7} + new_H{8} + new_H{9})/3;
%绘制频谱图
figure(3);
subplot(3,1,1);
plot(new_freq{1},r1_H);
title('r1频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 3000]);
subplot(3,1,2);
plot(new_freq{4},r4_H);
title('r4频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 3000]);
subplot(3,1,3);
plot(new_freq{7},r6_H);
title('r6频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 3000]);
%%
%检查r1_H,r4_H,r6_H中的非零元素个数
count1 = 0;
for i = 1:length(r1_H)
    if r1_H(i) ~= 0
        count1 = count1 + 1;
    end
end
count2 = 0;
for i = 1:length(r4_H)
    if r4_H(i) ~= 0
        count2 = count2 + 1;
    end
end
count3 = 0;
for i = 1:length(r6_H)
    if r6_H(i) ~= 0
        count3 = count3 + 1;
    end
end
%%
%求出三组r1_H,r4_H,r6_H非零元素对应的索引最大值
index1 = 0;
index2 = 0;
index3 = 0;
for i = 1:length(r1_H)
    if r1_H(i) ~= 0
        index1 = i;
    end
end
for i = 1:length(r4_H)
    if r4_H(i) ~= 0
        index2 = i;
    end
end
for i = 1:length(r6_H)
    if r6_H(i) ~= 0
        index3 = i;
    end
end
%打印索引处的频率
fprintf("r1_H的最大值对应的频率为：%f\n",new_freq{1}(index1));
fprintf("r4_H的最大值对应的频率为：%f\n",new_freq{4}(index2));
fprintf("r6_H的最大值对应的频率为：%f\n",new_freq{7}(index3));
%打印索引处的r1_H,r4_H,r6_H的值
fprintf("r1_H的最大值为：%f\n",r1_H(index1));
fprintf("r4_H的最大值为：%f\n",r4_H(index2));
fprintf("r6_H的最大值为：%f\n",r6_H(index3));

%截至频率定义为三个索引最大的一个的频率
index_max = max([index1,index2,index3]);
cut_freq = max([new_freq{1}(index1),new_freq{4}(index2),new_freq{7}(index3)]);
fprintf("截至频率为：%f\n",cut_freq);
%%
%根据截至频率设计低通滤波器
%截至频率为cut_freq
%采样频率为fs
%巴特沃斯低通滤波器
%阶数为4
%通带最大衰减为1dB
%通带截至频率为cut_freq
%阻带最小衰减为60dB

% %归一化截至角频率
% cut_freq = cut_freq/fs*pi;
% %设计低通滤波器
% [b,a] = butter(4,cut_freq,'low');
% %绘制滤波器的频率响应
% figure(4);
% freqz(b,a);
%%
%把三个频率截断在cut_freq处
%把r1_H,r4_H,r6_H截断在cut_freq处
%把new_freq{1},new_freq{4},new_freq{7}截断在cut_freq处
r1_H = r1_H(1:index_max);
r4_H = r4_H(1:index_max);
r6_H = r6_H(1:index_max);
new_freq{1} = new_freq{1}(1:index_max);
new_freq{4} = new_freq{4}(1:index_max);
new_freq{7} = new_freq{7}(1:index_max);
%%绘制频谱图
figure(5);
subplot(3,1,1);
plot(new_freq{1},r1_H);
title('r1频谱图');
xlabel('频率');
ylabel('幅度');
%横坐标调整到设定的截至频率
xlim([0 cut_freq]);
subplot(3,1,2);
plot(new_freq{4},r4_H);
title('r4频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 cut_freq]);
subplot(3,1,3);
plot(new_freq{7},r6_H);
title('r6频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 cut_freq]);
%%
%按index的步长为10扫描序列，只保留最大值，其余置零
for i = 1:20:length(r1_H)
    max1 = 0;
    for j = i:min(i+19,length(r1_H))
        if r1_H(j) > max1
            max1 = r1_H(j);
        end
    end
    for j = i:min(i+19,length(r1_H))
        if r1_H(j) ~= max1
            r1_H(j) = 0;
        end
    end
end
for i = 1:20:length(r4_H)
    max1 = 0;
    for j =  i:min(i+19,length(r1_H))
        if r4_H(j) > max1
            max1 = r4_H(j);
        end
    end
    for j =  i:min(i+19,length(r1_H))
        if r4_H(j) ~= max1
            r4_H(j) = 0;
        end
    end
end
for i = 1:20:length(r6_H)
    max1 = 0;
    for j = i:min(i+19,length(r1_H))
        if r6_H(j) > max1
            max1 = r6_H(j);
        end
    end
    for j =  i:min(i+19,length(r1_H))
        if r6_H(j) ~= max1
            r6_H(j) = 0;
        end
    end
end
%绘制频谱图
figure(6);
subplot(3,1,1);
plot(new_freq{1},r1_H);
title('r1频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 cut_freq]);
subplot(3,1,2);
plot(new_freq{4},r4_H);
title('r4频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 cut_freq]);
subplot(3,1,3);
plot(new_freq{7},r6_H);
title('r6频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0 cut_freq]);

%统计·r1_H,r4_H,r6_H中的非零元素个数
count1 = 0;
for i = 1:length(r1_H)
    if r1_H(i) ~= 0
        count1 = count1 + 1;
    end
end
count2 = 0;
for i = 1:length(r4_H)
    if r4_H(i) ~= 0
        count2 = count2 + 1;
    end
end
count3 = 0;
for i = 1:length(r6_H)
    if r6_H(i) ~= 0
        count3 = count3 + 1;
    end
end
%打印非零元素个数
fprintf("r1_H的非零元素个数为：%d\n",count1);
fprintf("r4_H的非零元素个数为：%d\n",count2);
fprintf("r6_H的非零元素个数为：%d\n",count3);
%%
%定义kong1,kong4,kong6分别拼接r1_H,new_freq{1},r4_H,new_freq{4},r6_H,new_freq{7}中的非零元素
kong1 = [];
index1 = [];%保存r1_H中非零元素的索引
kong4 = [];
index4 = [];%保存r4_H中非零元素的索引
kong6 = [];
index6 = [];%保存r6_H中非零元素的索引
%求解索引
for i = 1:length(r1_H)
    if r1_H(i) ~= 0
        index1 = [index1,i];
    end
end
for i = 1:length(r4_H)
    if r4_H(i) ~= 0
        index4 = [index4,i];
    end
end
for i = 1:length(r6_H)
    if r6_H(i) ~= 0
        index6 = [index6,i];
    end
end
%把索引指向的元素拼接到kong1,kong4,kong6中
%先接入new_freq{1}
for i = 1:length(index1)
    kong1 = [kong1,new_freq{1}(index1(i))];
end
for i = 1:length(index4)
    kong4 = [kong4,new_freq{4}(index4(i))];
end
for i = 1:length(index6)
    kong6 = [kong6,new_freq{7}(index6(i))];
end
%再接入r1_H,r4_H,r6_H
for i = 1:length(index1)
    kong1 = [kong1,r1_H(index1(i))];
end
for i = 1:length(index4)
    kong4 = [kong4,r4_H(index4(i))];
end
for i = 1:length(index6)
    kong6 = [kong6,r6_H(index6(i))];
end
%%
%把kong1,kong4,kong6保存到share.mat中
save('share.mat','kong1','kong4','kong6');

%%
function new_signal = processSignal(signal, freqs)
    % 初始化新信号
    new_signal = signal;
    
    % 计算误差范围内的数据点数
    n = round(50/(freqs(2) - freqs(1)));
    %对里面的部分取整数

    % 用于加快均值计算速度的初始和
    sum = 0;
    for i = 1:n
        sum = sum + signal(i);
    end
    
    % 凸峰幅值需为周围数据点均值的3倍以上
    mul = 3;
    
    % 存储上一个预选凸峰的索引值
    tmp = 1;
  
    % 遍历信号
    for i = 1:length(signal)
        sum = sum + signal(i);
  
        % 如果信号小于最大信号的1/15，则置零
        if signal(i) < max(signal) / 15
            new_signal(i) = 0;
            continue;
        end
  
        % 快速计算均值
        if i < n
            mean = (sum + signal(n + i)) / (n + i + 1);
        elseif i < length(signal) - n
            mean = (sum + signal(n + i) - signal(i - n - 1)) / (2 * n);
        else
            mean = (sum - signal(i - n - 1)) / (length(signal) - i + n);
        end
  
        % 检查凸峰
        if signal(i) > mean * mul
            if freqs(i) - freqs(tmp) < 20
                % 对间距小于20Hz的凸峰进行筛选
                if signal(tmp) < signal(i)
                    new_signal(tmp) = 0;
                    tmp = i;
                end
            end
        else
            new_signal(i) = 0;
        end
    end
  
    % 凸峰筛选
    threshold = max(new_signal) * 0.1;
    peaks = new_signal > threshold;
    new_signal = new_signal .* peaks;

    return 
end

%%
%截断函数
function [x1,x2] = cut(x)
    threshold1 = 1e-7;
    threshold2 = 0.08;
    pointcut = 0;
    point1 = 0;
    point2 = 0;
    point3 = 0;
    for i = 1:length(x)
        if abs(x(i)) < threshold1
            point1 = i;
            break;
        end
    end
    for i = point1:length(x)
        if abs(x(i)) > threshold2
            for j = i:length(x)
                if abs(x(j)) < threshold1
                    point2 = j;
                    break;
                end
            end
        end
    end
    for i = point2:length(x)
        if abs(x(i)) > threshold2
            for j = i:length(x)
                if abs(x(j)) < threshold1
                    point3 = j;
                    break;
                end
            end
        end
    end
    pointcut = point2;
    x1 = x(1:pointcut);
    x2 = x(pointcut+1:length(x));
end
