% 设置音频文件所在的文件夹路径
myFolder = 'E:\dsp大作业\data\面板敲击音'; % 替换为您音频文件夹的实际路径

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
savepath = 'E:\dsp大作业\picture\面板敲击音';
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
    saveas(gcf,fullfile(savepath,baseFileName(1:end-4)),'png');

    %截断Y1到index_max处
    Y1=Y1(1:index_max);
    %把截断的Y1保存到new_H中
    new_Y{k} = Y1;

end
%%
%对每个new_Y{i}进行三次独立的聚类，分别为count1，count2，count3类
%初始化一个二维数组，第一维表示文件序号，第二维表示聚类序号，数组元素表示对应的幅度值
new_Y11 = [];
new_Y12 = [];
new_Y13 = [];

%设置图片保存路径
savepath1 = 'E:\dsp大作业\picture\面板敲击音\聚类后的频谱图1';
savepath2 = 'E:\dsp大作业\picture\面板敲击音\聚类后的频谱图2';
savepath3 = 'E:\dsp大作业\picture\面板敲击音\聚类后的频谱图3';


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
    new_Y11{i} = zeros(1,2*count1);
    for j = 1:count1
        new_Y11{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y11{i}(j+count1) = C(j);
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
    new_Y12{i} = zeros(1,2*count2);
    for j = 1:count2
        new_Y12{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y12{i}(j+count2) = C(j);
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
    new_Y13{i} = zeros(1,2*count3);
    for j = 1:count3
        new_Y13{i}(j) = nearestPoints(j)*fs/length(new_Y{i});
        new_Y13{i}(j+count3) = C(j);
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
save('share.mat','new_Y11','new_Y12','new_Y13','-append');
%%
%把count1，count2，count3和cut_freq添加到share1.mat中
save('share1.mat','count1','count2','count3','cut_freq');
