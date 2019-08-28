cd('/mnt/f2f2edcb-b38e-4be6-869e-ce30964b605f/thesis_data/smt-hs-4-9')
ls = dir('*.csv');
av_im = zeros(512,513);
for f=1:length(ls)
    filename = ls(f).name;
    im = csvread(filename);
    av_im = im + av_im;
end    
%%
figure
    imagesc([0 3.05],[],im_10);
    set(gca,'YDir','normal')
    ax = gca;
    ax.YTickMode = 'manual';
    ax.YTick = [65   150   214   278   363   428   492];
    ax.YTickLabel = {'0.02','0.05','0.1','0.2','0.5','1','2'};
    hold on
    colorbar
    t=[];
    for i=0.01:0.01:2.5

        t=[t log_convert(i)];

    end
    i=0.01:0.01:2.5;
    plot(i,t,'r');
    xlabel('Logitudinal coefficient [$\mu m^2$/ms]','Interpreter','Latex')
    ylabel('Transverse coefficient [$\mu m^2$/ms]','Interpreter','Latex')