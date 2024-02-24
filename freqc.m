% 设置音频文件所在的文件夹路径
myFolder1 = 'E:\dsp大作业\data\空弦音'; % 替换为您音频文件夹的实际路径
myFolder2 = 'E:\dsp大作业\data\底板敲击音'; % 替换为您音频文件夹的实际路径


% 获取文件夹中所有.wav文件的列表
filePattern1 = fullfile(myFolder1, '*.wav'); % 替换为您的特定文件夹
wavFiles1 = dir(filePattern1);
filePattern2 = fullfile(myFolder2, '*.wav'); % 替换为您的特定文件夹
wavFiles2 = dir(filePattern2);

%读取空弦音和底板敲击音各自的第一个音频文件
baseFileName1 = wavFiles1(1).name;
fullFileName1 = fullfile(myFolder1, baseFileName1);
baseFileName2 = wavFiles2(1).name;
fullFileName2 = fullfile(myFolder2, baseFileName2);
[y1, Fs1] = audioread(fullFileName1);
[y2, Fs2] = audioread(fullFileName2);

%做FFT
Y1 = fft(y1);
Y2 = fft(y2);
%取前一半且只保留复数的模
Y1 = abs(Y1(1:length(Y1)/2));
Y2 = abs(Y2(1:length(Y2)/2));

%计算平均功率谱
P1 = abs(Y1).^2/length(Y1);
P2 = abs(Y2).^2/length(Y2);
P1_avg = mean(P1);
P2_avg = mean(P2);

%计算均方根值
rms1 = sqrt(mean(y1.^2));
rms2 = sqrt(mean(y2.^2));

p1_sum = sum(P1);
p2_sum = sum(P2);

%计算频谱质心
f1 = 0;
for i = 1:length(P1)
    f1 = f1 + i*P1(i);
end 
f1 = f1/p1_sum
f2 = 0;
for i = 1:length(P2)
    f2 = f2 + i*P2(i);
end
f2 = f2/p2_sum

%打印全部结果
fprintf('空弦音平均功率谱为：%f\n',P1_avg);
fprintf('底板敲击音平均功率谱为：%f\n',P2_avg);
fprintf('空弦音的均方根值为：%f\n',rms1);
fprintf('底板敲击音的均方根值为：%f\n',rms2);
fprintf('空弦音的频谱质心为：%f\n',f1);
fprintf('底板敲击音的频谱质心为：%f\n',f2);

