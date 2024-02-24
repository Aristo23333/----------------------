%遍历new_Y21,new_Y22,new_Y23中的每一个向量，分别计算其与kong1,kong2,kong3的相似性
%定义相似性结果向量cos_kong1,cos_kong2,cos_kong3
cos_kong1di = [];
cos_kong4di = [];
cos_kong6di = [];
cos_kong1mian = [];
cos_kong4mian = [];
cos_kong6mian = [];

%计算new_Y21中每一个元素与kong1的相似性
for i = 1:1:length(new_Y21)
    cos_kong1di(i) = cosine_similarity(new_Y21{i},kong1);
end

%计算new_Y22中每一个元素与kong2的相似性
for i = 1:1:length(new_Y22)
    cos_kong4di(i) = cosine_similarity(new_Y22{i},kong4);
end

%计算new_Y23中每一个元素与kong3的相似性
for i = 1:1:length(new_Y23)
    cos_kong6di(i) = cosine_similarity(new_Y23{i},kong6);
end

%计算new_Y11中每一个元素与kong1的相似性
for i = 1:1:length(new_Y11)
    cos_kong1mian(i) = cosine_similarity(new_Y11{i},kong1);
end

%计算new_Y12中每一个元素与kong2的相似性
for i = 1:1:length(new_Y12)
    cos_kong4mian(i) = cosine_similarity(new_Y12{i},kong4);
end

%计算new_Y13中每一个元素与kong3的相似性
for i = 1:1:length(new_Y13)
    cos_kong6mian(i) = cosine_similarity(new_Y13{i},kong6);
end
%%
%把kong1计算的两个cos_kong1di,cos_kong1mian绘制成热力图的形式，颜色越深，相似性越高
%注意绘制成两条热力图，一条是cos_kong1di,一条是cos_kong1mian
%上下排版，上面是cos_kong1di,下面是cos_kong1mian
%绘制cos_kong1di
figure(1)
imagesc(cos_kong1di);
colormap(jet);
colorbar;
title('cos_kong1底板');

%绘制cos_kong1mian
figure(2)
imagesc(cos_kong1mian);
colormap(jet);
colorbar;
title('cos_kong1面板');
%%
%把kong2计算的两个cos_kong2di,cos_kong2mian绘制成热力图的形式，颜色越深，相似性越高
%注意绘制成两条热力图，一条是cos_kong2di,一条是cos_kong2mian
%上下排版，上面是cos_kong2di,下面是cos_kong2mian
%绘制cos_kong2di
figure(3)
imagesc(cos_kong4di);
colormap(jet);
colorbar;
title('cos_kong4底板');

%绘制cos_kong2mian
figure(4)
imagesc(cos_kong4mian);
colormap(jet);
colorbar;
title('cos_kong4面板');
%%
%把kong3计算的两个cos_kong3di,cos_kong3mian绘制成热力图的形式，颜色越深，相似性越高
figure(5)
imagesc(cos_kong6di);
colormap(jet);
colorbar;
title('cos_kong6底板');

%绘制cos_kong3mian
figure(6)
imagesc(cos_kong6mian);
colormap(jet);
colorbar;
title('cos_kong6面板');


%%
%定义一个数组保存六个结果中的最大值
max_cos = [];
%定义一个数组保存六个结果中的最大值的索引
max_index = [];
%定义一个数组储存平均值
average_cos = [];

sum_cos = 0;
max_cos_i = 0;
max_index_i = 0;

for i = 1:1:length(cos_kong1di)
    if cos_kong1di(i)> max_cos_i
        max_cos_i = cos_kong1di(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong1di(i);
end

max_cos(1) = max_cos_i;
max_index(1) = max_index_i;
average_cos(1) = sum_cos/length(cos_kong1di);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;

for i = 1:1:length(cos_kong4di)
    if cos_kong4di(i)> max_cos_i
        max_cos_i = cos_kong4di(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong4di(i);
end

max_cos(2) = max_cos_i;
max_index(2) = max_index_i;
average_cos(2) = sum_cos/length(cos_kong4di);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;

for i = 1:1:length(cos_kong6di)
    if cos_kong6di(i)> max_cos_i
        max_cos_i = cos_kong6di(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong6di(i);
end

max_cos(3) = max_cos_i;
max_index(3) = max_index_i;
average_cos(3) = sum_cos/length(cos_kong6di);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;

for i = 1:1:length(cos_kong1mian)
    if cos_kong1mian(i)> max_cos_i
        max_cos_i = cos_kong1mian(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong1mian(i);
end

max_cos(4) = max_cos_i;
max_index(4) = max_index_i;
average_cos(4) = sum_cos/length(cos_kong1mian);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;

for i = 1:1:length(cos_kong4mian)
    if cos_kong4mian(i)> max_cos_i
        max_cos_i = cos_kong4mian(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong4mian(i);
end

max_cos(5) = max_cos_i;
max_index(5) = max_index_i;
average_cos(5) = sum_cos/length(cos_kong4mian);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;

for i = 1:1:length(cos_kong6mian)
    if cos_kong6mian(i)> max_cos_i
        max_cos_i = cos_kong6mian(i);
        max_index_i = i;
    end
    sum_cos = sum_cos + cos_kong6mian(i);
end

max_cos(6) = max_cos_i;
max_index(6) = max_index_i;
average_cos(6) = sum_cos/length(cos_kong6mian);
max_cos_i = 0;
max_index_i = 0;
sum_cos = 0;






%%
%定义余弦相似性计算函数，输入两个向量，输出相似性
function [similarity] = cosine_similarity(vector1,vector2)
    similarity = dot(vector1,vector2)/(norm(vector1)*norm(vector2));
end