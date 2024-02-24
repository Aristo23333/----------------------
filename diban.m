% 设置音频文件所在的文件夹路径
myFolder = 'E:\dsp大作业\data\底板敲击音'; % 替换为您音频文件夹的实际路径

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
savepath = 'E:\dsp大作业\picture\底板敲击音';


%%
%归一化截至角频率
cut_freq = cut_freq/fs*pi;
%设计低通滤波器
[b,a] = butter(4,cut_freq,'low');
%绘制滤波器的频率响应
figure(1);
freqz(b,a);
title('滤波器的频率响应');

%%
%循环读取文件夹中的所有.wav文件，仿照上面的代码计算频谱图，并画出处理前后的频谱图，且把全部处理后的频率特性保存到一个数组中
%这里的数组是一个二维数组，第一维表示文件序号，第二维表示频率序号，数组元素表示对应的幅度值
new_freq = [];
new_Y = [];
for k = 1:length(wavFiles)
    baseFileName = wavFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    %读取音频文件
    [y,fs] = audioread(fullFileName);

    %做傅里叶变换，注意这里的y是复数，所以要取绝对值并且只取前一半
    Y = abs(fft(y));
    Y = Y(1:length(Y)/2);
    %滤波  
    y1 = filter(b,a,y);
    Y1 = abs(fft(y1));
    Y1 = Y1(1:length(Y1)/2);
    %画出滤波后的频谱图和时域图
    figure(2);
    subplot(2,1,1);
    %绘制滤波后的时域图
    plot(y1);
    title('滤波后的时域图');
    xlabel('时间');
    ylabel('幅度');
   
    subplot(2,1,2);
    plot(Y1);
    title('滤波后的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
%     %保存图片
%     saveas(gcf,fullfile(savepath,baseFileName(1:end-4)),'png');

    %截断Y1到index_max处
    Y1=Y1(1:index_max);
    %把截断的Y1保存到new_H中
    new_Y{k} = Y1;

end

%%
% 利用kmeans聚类new_Y{2}为34类

k = 34;
%转置new_Y{2}
new_Y{2} = new_Y{2}';
%利用kmeans函数聚类
[idx,C] = kmeans(new_Y{2},k);
%画出聚类后的中心点，画在一张图上
%找到k个中心点的在频率上的索引
new_Y{2} = new_Y{2}';
% 初始化一个向量来存储每个质心最近点的索引
nearestPoints = zeros(1, k);

% 对于每个质心，找到最近的点
for i = 1:k
    % 遍历所有点，找到值最接近质心的点
    minDistance = inf;
    for j = 1:length(new_Y{2})
        distance = abs(new_Y{2}(j) - C(i));
        if distance < minDistance
            minDistance = distance;
            nearestPoints(i) = j;
        end
    end
end

%重新绘制频谱图，在index_Y2处的值为中心点的值，其余为0
Y2 = zeros(1,length(new_Y{2}));
for i = 1:k
    Y2(nearestPoints(i)) = C(i);
end

%求出index_Y2对应的频率
freq_Y2 = nearestPoints*fs/length(new_Y{2});


%画出聚类后的频谱图
figure(4);
plot(Y2);
title('聚类后的频谱图');
xlabel('频率');
ylabel('幅度');
xlim([0,20000]);
%%
%对每个new_Y{i}进行三次独立的聚类，分别为count1，count2，count3类
%初始化一个二维数组，第一维表示文件序号，第二维表示聚类序号，数组元素表示对应的幅度值
new_Y21 = [];
new_Y22 = [];
new_Y23 = [];

%设置图片保存路径
savepath1 = 'E:\dsp大作业\picture\底板敲击音\聚类后的频谱图1';
savepath2 = 'E:\dsp大作业\picture\底板敲击音\聚类后的频谱图2';
savepath3 = 'E:\dsp大作业\picture\底板敲击音\聚类后的频谱图3';


for i = 1:length(new_Y)
    %转置new_Y{i}
    new_Y{i} = new_Y{i}';
    %利用kmeans函数聚类
    [idx,C] = kmeans(new_Y{i},count1);
    %画出聚类后的中心点，画在一张图上
    %找到count1个中心点的在频率上的索引
    new_Y{i} = new_Y{i}';
    % 初始化一个向量来存储每个质心最近点的索引
    nearestPoints = zeros(1, count1);

    % 对于每个质心，找到最近的点
    for j = 1:count1
        % 遍历所有点，找到值最接近质心的点
        minDistance = inf;
        for k = 1:length(new_Y{i})
            distance = abs(new_Y{i}(k) - C(j));
            if distance < minDistance
                minDistance = distance;
                nearestPoints(j) = k;
            end
        end
    end

    %重新绘制频谱图，在index_Y2处的值为中心点的值，其余为0
    Y2 = zeros(1,length(new_Y{i}));
    for j = 1:count1
        Y2(nearestPoints(j)) = C(j);
    end
    
    %取出index_Y2对应的频率和幅度，保存到new_Y21中其中new_Y21{i}的前count1个元素为频率，后count1个元素为幅度
    new_Y21{i} = zeros(1,2*count1);
    for j = 1:count1
        new_Y21{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y21{i}(j+count1) = C(j);
    end
    
    %画出聚类前后的频谱图
    figure(5);
    subplot(2,1,1);
    plot(new_Y{i});
    title('聚类前的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    subplot(2,1,2);
    plot(Y2);
    title('聚类后的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    %保存图片
    saveas(gcf,fullfile(savepath1,wavFiles(i).name(1:end-4)),'png');

end

for i = 1:length(new_Y)
    %转置new_Y{i}
    new_Y{i} = new_Y{i}';
    %利用kmeans函数聚类
    [idx,C] = kmeans(new_Y{i},count2);
    %画出聚类后的中心点，画在一张图上
    %找到count2个中心点的在频率上的索引
    new_Y{i} = new_Y{i}';
    % 初始化一个向量来存储每个质心最近点的索引
    nearestPoints = zeros(1, count2);

    % 对于每个质心，找到最近的点
    for j = 1:count2
        % 遍历所有点，找到值最接近质心的点
        minDistance = inf;
        for k = 1:length(new_Y{i})
            distance = abs(new_Y{i}(k) - C(j));
            if distance < minDistance
                minDistance = distance;
                nearestPoints(j) = k;
            end
        end
    end

    %重新绘制频谱图，在index_Y2处的值为中心点的值，其余为0
    Y2 = zeros(1,length(new_Y{i}));
    for j = 1:count2
        Y2(nearestPoints(j)) = C(j);
    end
    
    %取出index_Y2对应的频率和幅度，保存到new_Y22中其中new_Y22{i}的前count2个元素为频率，后count2个元素为幅度
    new_Y22{i} = zeros(1,2*count2);
    for j = 1:count2
        new_Y22{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y22{i}(j+count2) = C(j);
    end
    
    %画出聚类前后的频谱图
    figure(6);
    subplot(2,1,1);
    plot(new_Y{i});
    title('聚类前的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    subplot(2,1,2);
    plot(Y2);
    title('聚类后的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    %保存图片
    saveas(gcf,fullfile(savepath2,wavFiles(i).name(1:end-4)),'png');

end

for i = 1:length(new_Y)
    %转置new_Y{i}
    new_Y{i} = new_Y{i}';
    %利用kmeans函数聚类
    [idx,C] = kmeans(new_Y{i},count3);
    %画出聚类后的中心点，画在一张图上
    %找到count3个中心点的在频率上的索引
    new_Y{i} = new_Y{i}';
    % 初始化一个向量来存储每个质心最近点的索引
    nearestPoints = zeros(1, count3);

    % 对于每个质心，找到最近的点
    for j = 1:count3
        % 遍历所有点，找到值最接近质心的点
        minDistance = inf;
        for k = 1:length(new_Y{i})
            distance = abs(new_Y{i}(k) - C(j));
            if distance < minDistance
                minDistance = distance;
                nearestPoints(j) = k;
            end
        end
    end

    %重新绘制频谱图，在index_Y2处的值为中心点的值，其余为0
    Y2 = zeros(1,length(new_Y{i}));
    for j = 1:count3
        Y2(nearestPoints(j)) = C(j);
    end
    
    %取出index_Y2对应的频率和幅度，保存到new_Y23中其中new_Y23{i}的前count3个元素为频率，后count3个元素为幅度
    new_Y23{i} = zeros(1,2*count3);
    for j = 1:count3
        new_Y23{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y23{i}(j+count3) = C(j);
    end
    
    %画出聚类前后的频谱图
    figure(7);
    subplot(2,1,1);
    plot(new_Y{i});
    title('聚类前的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    subplot(2,1,2);
    plot(Y2);
    title('聚类后的频谱图');
    xlabel('频率');
    ylabel('幅度');
    xlim([0,20000]);
    %保存图片
    saveas(gcf,fullfile(savepath3,wavFiles(i).name(1:end-4)),'png');

end
%%
%把new_Y21，new_Y22，new_Y23添加到share.mat中
save('share.mat','new_Y21','new_Y22','new_Y23','-append');

%%
%定义x1和x2
x1 = 0;
x2 = 0;
%调用cut函数，截断b3.wav
[x1,x2] = cut(x);


%绘制截断后的频域图，注意只画出前半部分
X1 = abs(fft(x1));
X2 = abs(fft(x2));

figure(3);
subplot(2,1,1);
plot(X1(1:floor(length(X1)/2)));
title('截断后的频域图1');
xlabel('频率');
ylabel('幅度');
xlim([0,5000]);
subplot(2,1,2);
plot(X2(1:floor(length(X2)/2)));
title('截断后的频域图2');
xlabel('频率');
ylabel('幅度');
xlim([0,5000]);
%%
%定义x1_max和x2_max为x1和x2的最大值，并保存索引
x1_max = max(x1);
x2_max = max(x2);
%定义x1_max_index和x2_max_index为x1和x2的最大值的索引
x1_max_index = 0;
x2_max_index = 0;
%找到x1和x2的最大值的索引
for i = 1:length(x1)
    if x1(i) == x1_max
        x1_max_index = i;
        break;
    end
end
for i = 1:length(x2)
    if x2(i) == x2_max
        x2_max_index = i;
        break;
    end
end
%最大值以外的部分置零
x1_input = zeros(1,length(x1));
x2_input = zeros(1,length(x2));

%最大值索引处的值置为最大值
x1_input(x1_max_index) = x1_max;
x2_input(x2_max_index) = x2_max;

%画出x1_input1和x1_input2的时域图
figure(4);
subplot(2,1,1);
plot(x1_input);
title('x1_input的时域图');
xlabel('时间');
ylabel('幅度');
subplot(2,1,2);
plot(x2_input);
title('x2_input的时域图');
xlabel('时间');
ylabel('幅度');
%%
%画出x1_input1和x1_input2的频域图
X1_input = abs(fft(x1_input));
X2_input = abs(fft(x2_input));

figure(5);
subplot(2,1,1);
plot(X1_input(1:floor(length(X1_input)/2)));
title('x1_input的频域图');
xlabel('频率');
ylabel('幅度');
subplot(2,1,2);
plot(X2_input(1:floor(length(X2_input)/2)));
title('x2_input的频域图');
xlabel('频率');
ylabel('幅度');

%%
%定义x1_signal和x2_signal
x1_signal = X1 ;
x2_signal = X2 ;
%定义x1_freqs和x2_freqs,表示频率
x1_freqs = (0:length(x1_signal)-1)*fs/length(x1_signal);
x2_freqs = (0:length(x2_signal)-1)*fs/length(x2_signal);
%处理信号
x1_signal = processSignal(x1_signal, x1_freqs);
x2_signal = processSignal(x2_signal, x2_freqs);
%画出x1_signal和x2_signal的频域图
figure(6);
subplot(2,1,1);
plot(x1_signal);
title('x1_signal的频域图');
xlabel('频率');
ylabel('幅度');
xlim([0,5000]);
subplot(2,1,2);
plot(x2_signal);
title('x2_signal的频域图');
xlabel('频率');
ylabel('幅度');
xlim([0,5000]);

%%
%把上面的部分封装成函数，输入x，输出x1和x2
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

%%

