
% load('data/contrast_result4.mat')
colorsTab = [
    0.0, 0.0, 0.0;     % 黑色
    0.49, 0.18, 0.56;  % 紫色
    1.0, 0.0, 0.0;     % 醒目的红色
    1.0, 0.0, 1.0;     % 醒目的洋红色
    0.0, 1.0, 1.0;     % 醒目的浅蓝色
    0.0, 1.0, 0.0;     % 醒目的绿色
];
filtername=["TEKF", "NEKF", "PICDAKF", "PSCDAKF", "ICDAKF", "SCDAKF"];
t = 1:ts;
linewidth_prmse = 2;
linewidth_vrmse = 2;
fontsize = 16;

figure(1)
subplot(1,2,1)
for i = 1:length(filtername)
    hold on
    h = plot(t, smooth(pos_RMSE(:,i),50),'color', colorsTab(i,:), 'LineWidth', linewidth_prmse);
    L(i) = legend(h, filtername(i), 'FontSize', 16);
    legend off;
end
legend show;
% axis([0,10000,3,6.6])
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
xlabel('Time(s)','FontSize',24, 'FontName', 'Times New Roman', 'FontWeight', 'bold')
ylabel('PRMSE(m)','FontSize',24, 'FontName', 'Times New Roman', 'FontWeight', 'bold')

subplot(1,2,2)
for i = 1:length(filtername)
    hold on
    h = plot(t, smooth(vel_RMSE(:,i),50),'color', colorsTab(i,:), 'LineWidth', linewidth_vrmse);
    L(i) = legend(h, filtername(i));
    legend off;
end
legend show;
set(gca, 'FontSize', 20, 'FontWeight', 'bold')
xlabel('Time(s)','FontSize',24, 'FontName', 'Times New Roman', 'FontWeight', 'bold')
ylabel('VRMSE(m)','FontSize',24, 'FontName', 'Times New Roman', 'FontWeight', 'bold')

figure(2)
for i = 1:length(filtername)-2
    hold on
    h = plot(t, smooth(saveale(:,i)),'color', colorsTab(i+2,:), 'LineWidth', 2);
    L(i) = legend(h, filtername(i+2));
    legend off;
end
legend show;
set(gca, 'FontSize', 24, 'FontWeight', 'bold')
xlabel('Time(s)','FontSize',28, 'FontName', 'Times New Roman', 'FontWeight', 'bold')
ylabel('ALE','FontSize',28, 'FontName', 'Times New Roman', 'FontWeight', 'bold')



