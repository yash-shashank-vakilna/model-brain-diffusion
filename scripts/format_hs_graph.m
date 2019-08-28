function format_hs_graph(full_file_path)
im = csvread(full_file_path);
%f = figure('visible','off');
f = figure;
imagesc([0 3.05],[],im);
set(gca,'YDir','normal')
ax = gca;
ax.YTickMode = 'manual';
ax.YTick = [65   150   214   278   363   428   492];
ax.YTickLabel = {'0.02','0.05','0.1','0.2','0.5','1','2'};
hold on
t = zeros(1,250); n=1;
for i=0.01:0.01:2.5
    
    t(n)= log_convert(i);
    n=n+1;
end
i=0.01:0.01:2.5;
plot(i,t,'r');
xlabel('logitudinal coefficient')
ylabel('transverse coefficient')
colorbar
full_file_path=replace(full_file_path,'.csv','');
saveas(f,full_file_path,'png');
disp(full_file_path) 
end
