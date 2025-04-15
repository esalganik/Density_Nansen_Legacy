% data processing
close all; clc; clear;
project = "NL_density.xlsx";
T1 = readtable(project,'Sheet','density');
T2 = readtable(project,'Sheet','temperature');

n = table2array(T1(:,1)); rho_all = table2array(T1(:,19)); t_all = table2array(T1(:,2)); zzrho_all = table2array(T1(:,9:10)); T_lab_all = table2array(T1(:,11)); Srho_all = table2array(T1(:,12));
hrho_all = table2array(T1(:,8)); drho_all = table2array(T1(:,8))-table2array(T1(:,6)); vg_orig_all = table2array(T1(:,14)); lon_all = table2array(T1(:,4)); lat_all = table2array(T1(:,5));  
n2 = table2array(T2(:,1)); hT_all = table2array(T2(:,8)); dT_all = table2array(T2(:,8))-table2array(T2(:,6)); zT_all = table2array(T2(:,9)); T_all = table2array(T2(:,11)); 

for i = 1:max(n)
    rho{i} = rho_all(n == i); zzrho{i} = zzrho_all(n == i,:); Srho{i} = Srho_all(n == i); T_lab{i} = T_lab_all(n == i)*0-20; t0{i} = t_all(n == i); t(i) = t0{i}(1);
    drho{i} = drho_all(n == i); vg_orig{i} = vg_orig_all(n == i); hrho0{i} = hrho_all(n == i); hrho(i) = hrho0{i}(1);
    lon0{i} = lon_all(n == i); lat0{i} = lat_all(n == i); lon(i) = lon0{i}(1); lat(i) = lat0{i}(1);
    T{i} = T_all(n2 == i); zT{i} = zT_all(n2 == i); dT{i} = dT_all(n2 == i); hT0{i} = hT_all(n2 == i); hT(i) = hT0{i}(1); 
end
for i = 31:34; t0{i} = t0{i}-365; t(i) = t(i)-365; end
clearvars -except n hT hrho dT drho zT zzrho T T_lab Srho rho t t0 hrho0 lon lat c

for i = 1:max(n)
    zrho{i} = mean(zzrho{i},2) * hrho(i) / zzrho{i}(end,2);
    T_rho{i} = interp1(zT{i} * hT(i) / zT{i}(end),T{i},zrho{i},'linear','extrap');
    F1_pr_rho = -4.732-22.45*T_lab{i} - 0.6397*T_lab{i}.^2 - 0.01074*T_lab{i}.^3;
    F2_pr_rho = 8.903*10^-2 - 1.763*10^-2*T_lab{i} - 5.33*10^-4*T_lab{i}.^2 - 8.801*10^-6*T_lab{i}.^3;
    vb_pr_rho = rho{i} .* Srho{i} ./ F1_pr_rho; % brine volume for T_lab
    rhoi_pr = (917-1.403*10^-1*T_lab{i}); % pure ice density, Pounder (1965)
    vg_pr{i} = max((1-rho{i}.*(F1_pr_rho-rhoi_pr.*Srho{i}/1000.*F2_pr_rho)./(rhoi_pr.*F1_pr_rho)),0,'includenan'); % gas volume for T_lab
    F3_pr = rhoi_pr.*Srho{i}/1000./(F1_pr_rho-rhoi_pr.*Srho{i}/1000.*F2_pr_rho);
    rhoi_rho = (917-1.403*10^-1*T_rho{i}); % pure ice density, Pounder (1965) for T_insitu
    F1_rho = -4.732-22.45*T_rho{i} - 0.6397*T_rho{i}.^2 - 0.01074*T_rho{i}.^3; % Cox and Weeks (1983)
    F1_rho(T_rho{i}>-2) = -4.1221*10^-2 + -1.8407*10^1*T_rho{i}(T_rho{i}>-2).^1 + 5.8402*10^-1*T_rho{i}(T_rho{i}>-2).^2 + 2.1454*10^-1*T_rho{i}(T_rho{i}>-2).^3; % F1 from Lepparanta and Manninen (1988)
    F2_rho = 8.903*10^-2 - 1.763*10^-2*T_rho{i} - 5.33*10^-4*T_rho{i}.^2 - 8.801*10^-6*T_rho{i}.^3; % Cox and Weeks (1983)
    F2_rho(T_rho{i}>-2) = 9.0312*10^-2 + -1.6111*10^-2*T_rho{i}(T_rho{i}>-2).^1 + 1.2291*10^-4*T_rho{i}(T_rho{i}>-2).^2 + 1.3603*10^-4*T_rho{i}(T_rho{i}>-2).^3; % F2 from Lepparanta and Manninen (1988)
    F3_rho = rhoi_rho.*Srho{i}/1000./(F1_rho-rhoi_rho.*Srho{i}/1000.*F2_rho);
    vb_rho{i} = vb_pr_rho .* F1_pr_rho ./ F1_rho / 1000; vb_rho{i}(vb_rho{i} > 0.6) = NaN; vb_rho{i}(vb_rho{i} < 0) = NaN;  % Brine volume for T_insitu, CW + LM
    vg{i} = max(0,(1-(1-vg_pr{i}).*(rhoi_rho./rhoi_pr).*(F3_pr.*F1_pr_rho./F3_rho./F1_rho)),'includenan'); % Gas volume for T_insitu, CW + LM
    rho_si{i} = (-vg_pr{i}+1).*rhoi_rho.*F1_rho./(F1_rho-rhoi_rho.*Srho{i}/1000.*F2_rho); rho_si{i}(isnan(vb_rho{i})) = NaN; rho_si{i}(isnan(vb_rho{i})) = NaN; % density
    rho_si_bulk(i) = nanmean(rho_si{i}); T_bulk(i) = nanmean(T_rho{i});
end
clearvars -except hT hrho dT drho zT zrho T T_lab Srho rho rho_si vb_rho vg vg_pr t t0 lon lat n rho_si_bulk T_rho T_bulk c

%% Nansen legacy overview: map, density vs time, density vs temperature (6 x 2 in)
figure
tile = tiledlayout(1,3); tile.TileSpacing = 'compact'; tile.Padding = 'none'; ax = gobjects(1,3);
ax(1) = nexttile; % bathymetry
m_proj('lambert','lons',[-20 40],'lat',[77 88]);
[cs,~]=m_etopo2('contourf',-6000:100:0,'edgecolor','none');
m_gshhs_l('patch',[.7 .7 .7],'edgecolor','none'); % coastline
m_grid('linewi',1,'layer','top','FontSize',8); % axis settings
clim([-6000 0]);
colormap(ax(1),m_colmap('blue'));
[ax,~]=m_contfbar(.65,[1.0 1.8],[-6000 0],-6000:100:0,'edgecolor','none','endpiece','no','axfrac',.03,'tickdir','out'); title(ax,'Depth (m)','fontsize',6,'fontweight','normal'); % vert. colorbar settings
set(ax,'fontsize',6);
load("lipari100.mat");
for i = 1:max(n)
    m_line(lon(i),lat(i),'marker','o','color',lipari100(round(hrho(i)/max(hrho),1)*100,:),'markerfacecolor',lipari100(round(hrho(i)/max(hrho),1)*100,:),'linewi',0.1,'markersize',3.5);
end
p = m_text(16,79,'Svalbard'); set(p,'HorizontalAlignment','center','FontSize',7);

ax(2) = nexttile;
for i = 1:max(n)
    p1 = plot(t0{i},rho_si{i},'o','Color',[.7 .7 .7]); p1.MarkerSize = 1.2; set(p1,'markerfacecolor',get(p1,'color')); hold on
end
load("lipari100.mat");
for i = 1:max(n)
    p2 = plot(t(i),rho_si_bulk(i),'o','Color',lipari100(round(hrho(i)/max(hrho),1)*100,:)); p2.MarkerSize = 3.0; set(p2,'markerfacecolor',get(p2,'color')); hold on
end
leg = legend([p1 p2],'Sections','Full thickness','box','off'); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.5,18*0.5];
hXLabel = xlabel('Time'); hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
colormap(ax(2),lipari100); clim([0 1.8]);
ax = gca; ax.XTick = datetime(['01-Feb-2021';'01-Mar-2021';'01-Apr-2021';'01-May-2021';'01-Jun-2021';'01-Jul-2021';'01-Aug-2021';'01-Sep-2021';'01-Oct-2021']); datetick('x','mmm','keepticks'); xtickangle(0); % time

ax(3) = nexttile;
for i = 1:max(n)
    p1 = plot(T_rho{i},rho_si{i},'o','Color',[.7 .7 .7]); p1.MarkerSize = 1.2; set(p1,'markerfacecolor',get(p1,'color')); hold on
end
for i = 1:max(n)
    p2 = plot(T_bulk(i),rho_si_bulk(i),'o','Color',lipari100(round(hrho(i)/1.8,1)*100,:)); p2.MarkerSize = 3.0; set(p2,'markerfacecolor',get(p2,'color')); hold on
end
hXLabel = xlabel('Ice temperature (°C)'); hYLabel = ylabel('Ice density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
colormap(ax(3),lipari100); clim([0 1.8]);
hBar1 = colorbar; ylabel(hBar1,'Ice thickness (m)','FontSize',7);
leg = legend([p1 p2],'Sections','Full thickness','box','off'); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.5,18*0.5];

annotation('textbox',[0 .51 0.02 .51],'String','(a)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.20 .51 0.21 .51],'String','(b)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
annotation('textbox',[0.40 .51 0.42 .51],'String','(c)','FontSize',8,'EdgeColor','none','HorizontalAlignment','center');
clearvars ax tile p p1 p2 leg lipari100 hXLabel hYLabel gca cs i 

%% Nansen Legacy vs historical observations
load("density_historical.mat"); % data import

figure
tile = tiledlayout(1,2); tile.TileSpacing = 'compact'; tile.Padding = 'none';
nexttile
p = plot(t_hist([4 6:32 34:end]),rho_hist([4 6:32 34:end]),'kx'); p.MarkerSize = 4; hold on % Timco & Frederking 1996, doi:10.1016/0165-232X(95)00007-X
p = plot(t_leg5,rho_leg5,'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, leg 5, updated
p = plot(t-365,rho_si_bulk,'o','Color',c{2},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % Nansen Legacy
p = plot(t_fyi_mosaic(1:9)+365,rho_fyi_mosaic(1:9),'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, updated
p = plot(t_fyi_mosaic(10:23),rho_fyi_mosaic(10:23),'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, updated
leg = legend('historical','MOSAiC','Nansen Legacy','','','box','off','NumColumns',1); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.4,18*0.4];
hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hYLabel gca],'FontSize',8,'FontWeight','normal'); ylim([860 930]);
ax = gca; ax.XTick = datetime(['01-Jan-2020';'01-Mar-2020';'01-May-2020';'01-Jul-2020';'01-Sep-2020';'01-Nov-2020';'01-Jan-2021']); datetick('x','mmm','keepticks'); xtickangle(0); % time

nexttile
p = plot(T_hist([1:32 34:end]),rho_hist([1:32 34:end]),'kx'); p.MarkerSize = 4; hold on % Timco and Frederking, 1996
p = plot(T_fyi_mosaic,rho_fyi_mosaic,'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, updated
p = plot(T_leg5,rho_leg5,'o','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC, updated
p = plot(T_bulk,rho_si_bulk,'o','Color',c{2},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % Nansen Legacy
leg = legend('historical','MOSAiC','','Nansen Legacy','box','off','NumColumns',1); set(leg,'FontSize',6,'Location','southwest'); leg.ItemTokenSize = [30*0.7,18*0.7];
hXLabel = xlabel('FYI temperature (°C)'); hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
ylim([860 930]);
annotation('textbox',[0 .51 0.02 .51],'String','(a)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.34 .51 0.35 .51],'String','(b)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
clearvars ax hXLabel hYLabel i leg p T_anal rho_si_cox tile

%% Relationship with ice thickness, parametrization vs ice temperature
figure
tile = tiledlayout(1,2); tile.TileSpacing = 'compact'; tile.Padding = 'none';
nexttile
load("density_historical.mat"); % data import
load("lipari100.mat");
for i = [5:6 8:9 11:23]
    p1 = plot(T_fyi_mosaic(i),rho_fyi_mosaic(i),'<','Color',lipari100(round(h_fyi_mosaic(i)/100/2.5,1)*100,:)); p1.MarkerSize = 3.5; set(p1,'markerfacecolor',get(p1,'color')); hold on
end
for i = 1:length(rho_leg5)
    p = plot(T_leg5(i),rho_leg5(i),'<','Color',lipari100(round(h_leg5(i)/2.5,1)*100,:)); p.MarkerSize = 3.5; set(p,'markerfacecolor',get(p,'color')); hold on
end
for i = 1:max(n)
    p2 = plot(T_bulk(i),rho_si_bulk(i),'o','Color',lipari100(round(hrho(i)/2.5,1)*100,:)); p2.MarkerSize = 3.5; set(p2,'markerfacecolor',get(p2,'color')); hold on
end
hXLabel = xlabel('Ice temperature (°C)'); hYLabel = ylabel('Ice density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
colormap(lipari100); clim([0 2.5]);
hBar1 = colorbar; ylabel(hBar1,'Ice thickness (m)','FontSize',8);
leg = legend([p1 p2],'MOSAiC','Nansen Legacy','box','off'); set(leg,'FontSize',8,'Location','southwest'); leg.ItemTokenSize = [30*0.5,18*0.5];
xlim([-10 0]); ylim([868 920]);

% Density parametrization
offset = 914.1;
x = [T_fyi_mosaic([5:6 8:9 11:13 15:23]) T_bulk];
y = [rho_fyi_mosaic([5:6 8:9 11:13 15:23]) rho_si_bulk]-offset;

[xData,yData] = prepareCurveData(x,y);
ft = fittype('exp1'); opts = fitoptions('Method','NonlinearLeastSquares'); opts.Display = 'Off'; opts.Normalize = 'on'; opts.Robust = 'LAR'; % Set up fittype and options.
[fitresult,gof] = fit(xData,yData,ft,opts);
Xfit=linspace(min(xData),max(xData),100); Yfit=fitresult(Xfit);

nexttile
p = plot(T_fyi_mosaic([5:6 8:9 11:13 15:23]),rho_fyi_mosaic([5:6 8:9 11:13 15:23]),'<','Color',c{1},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); hold on
p = plot(T_leg5,rho_leg5,'<','Color',c{3},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % MOSAiC
p = plot(T_bulk,rho_si_bulk,'o','Color',c{2},'LineWidth',2); p.MarkerSize = 2; set(p,'markerfacecolor',get(p,'color')); % Nansen Legacy
plot(Xfit,Yfit+offset,'Color',c{4},'LineWidth',2); % exponential fit
p = text(-9.5,907,'914.1 - 2.8 exp(1.8T)'); set(p,'Color',c{4},'HorizontalAlignment','left','FontSize',9);
leg = legend('MOSAiC, Nov-Jul','MOSAiC, Sep','Nansen Legacy','Exponential fit (R2 = 0.97)','box','off','NumColumns',1); set(leg,'FontSize',8,'Location','southwest'); leg.ItemTokenSize = [30*0.7,18*0.7];
hXLabel = xlabel('Ice temperature (°C)'); hYLabel = ylabel('FYI density (kg  m^-^3)'); set([hXLabel hYLabel gca],'FontSize',8,'FontWeight','normal');
xlim([-10 0]); ylim([868 920]);
annotation('textbox',[0 .51 0.02 .51],'String','(a)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
annotation('textbox',[0.38 .51 0.38 .51],'String','(b)','EdgeColor','none','HorizontalAlignment','center','FontSize',8);
