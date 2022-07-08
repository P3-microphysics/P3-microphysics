
clear all

rrr = 6000;
ccc = rrr-49;

% 


load ./outQcld.dat
t00 = outQcld(1:50,6);

t0 = outQcld(ccc:rrr,6);
% td0 = outQcld(ccc:rrr,7);
vmi0 = outQcld(ccc:rrr,12);
z0 = outQcld(ccc:rrr,3);
z20 = flipud(z0);
qs0 = outQcld(ccc:rrr,10)*10.^-8;
qg0 = outQcld(ccc:rrr,8)*10.^-8;
qr0 = outQcld(ccc:rrr,7)*10.^-8;
qil0 = outQcld(ccc:rrr,9)*10.^-8;
qv0 = outQcld(ccc:rrr,11);
rho0 = outQcld(ccc:rrr,13);

frim0 = qg0./(qs0-qil0);
fliq0 = qil0./(qs0);
load ./outNcld.dat
ns0 = outNcld(ccc:rrr,5);
nr0 = outNcld(ccc:rrr,4);
zi0 = outNcld(ccc:rrr,9);
lami0 = outNcld(ccc:rrr,11);
mui0 = outNcld(ccc:rrr,10);
dm0 = outNcld(ccc:rrr,7);
ze0 = outNcld(ccc:rrr,12);

% 

cd ..
cd 1D_ORIG_case1
load ./outQcld.dat
t00 = outQcld(1:50,6);
t = outQcld(ccc:rrr,6);
% td = outQcld(ccc:rrr,7);
% vmi = outQcld(ccc:rrr,13);
z = outQcld(ccc:rrr,3);
z0 = outQcld(ccc:rrr,3);

z2 = flipud(z);
qs = outQcld(ccc:rrr,9)*10.^-8;
qg = outQcld(ccc:rrr,8)*10.^-8;
qr = outQcld(ccc:rrr,7)*10.^-8;
qv = outQcld(ccc:rrr,10);
% qc = outQcld(ccc:rrr,8);
rho = outQcld(ccc:rrr,12);
vmi = outQcld(ccc:rrr,11);
frim = outQcld(ccc:rrr,13);
load ./outNcld.dat
ns = outNcld(ccc:rrr,5);
nr = outNcld(ccc:rrr,4);
zi = outNcld(ccc:rrr,9);
lami = outNcld(ccc:rrr,10);
mui = outNcld(ccc:rrr,11);
dm = outNcld(ccc:rrr,7);
ze = outNcld(ccc:rrr,12);

% 
% figure(1)
%       plot(t00,z0*10.^-3,'Color','r','linestyle',':','LineWidth',3.5)          
%       hold on
%     plot(t0,z0*10.^-3,'Color','k','LineWidth',1)
%         hold on
%             plot(t,z*10.^-3,'Color','k','linestyle','--','LineWidth',3.5)
%     ylabel('z [km]')
%     xlabel('T [°C]')
%     set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
%     set(gca, 'xlim',[-3 2],'Xtick',[-3:1:2])
%     title('(a)')
%         leg = legend('T_initial','P3_MOD','P3_ORIG','Location','SE')
%         set(leg,'Interpreter','none')
%             set(gca,'Fontsize',20,'FontName','Arial','fontweight','bold') 
%  line('XData',[-3 2],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%  line('XData',[0 0],'YData',[0 1],'Linewidth',2.5,'Color','k','LineStyle',':')
% set(gcf,'Position',[10 10 400 400])

figure(2)
% plot(qr0*1000,z0*10.^-3,'Color',[0 .7 0],'LineWidth',1) 
%     hold on
    plot(qs0*1000,z0*10.^-3,'Color','k','lineStyle','-','LineWidth',1)
%     hold on
%         plot((qil0)*1000,z0*10.^-3,'Color','r','lineStyle','-','LineWidth',1)
%             hold on
%     plot(qr*1000,z*10.^-3,'Color',[0 .7 0],'linestyle','--','LineWidth',3.5) 
    hold on
    plot(qs*1000,z*10.^-3,'Color','k','lineStyle','--','LineWidth',3.5)
    ylabel('z [km]')
    xlabel('q [g kg^{-1}]')
    title('(b)')
    legend('q_{rain}','q_{i,tot}','q_{i,liq}','Location','NW')
        set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
                set(gca, 'xlim',[0 0.3],'XTick',[0:0.1:0.3])
            set(gca,'Fontsize',20,'FontName','Arial','fontweight','bold') 
    line('XData',[0 0.3],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
set(gcf,'Position',[10 10 400 400])
% 
    figure(3)
plot((ns0),z0*10.^-3,'Color','k','LineWidth',1)
% hold on
% plot(nr0,z0*10.^-3,'Color',[0 .7 0],'LineWidth',1)
hold on
plot((ns),z*10.^-3,'Color','k','linestyle','--','LineWidth',3.5)
% hold on
% plot(nr,z*10.^-3,'Color',[0 .7 0],'linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
        xlabel('N [kg^{-1}]')
        title('(c)')
                set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold') 
    legend('N_{i,tot}','N_{rain}','Location','SE')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
%                         set(gca, 'Xlim',[0 5000],'XTick',[0:2500:5000])
%         line('XData',[0 5000],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
set(gcf,'Position',[10 10 400 400])

figure(4)
plot(frim0,z0*10.^-3,'Color','b','Linestyle','-','LineWidth',1)
hold on
plot(fliq0,z0*10.^-3,'Color','k','Linestyle','-','LineWidth',1)
hold on
plot(frim,z*10.^-3,'Color','b','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    xlabel('Fraction')
    title('(d)')
                    set(gca,'Fontsize',20,'FontName','Arial','fontweight','bold') 
        set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    line('XData',[0 1],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
                set(gca, 'Xlim',[0 1],'XTick',[0:0.2:1])
set(gcf,'Position',[10 10 400 400])

                figure(5)
plot(rho0,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(rho,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
                        set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('\rho [kg m^{-3}]')
    title('(e)')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold') 
        line('XData',[0 1000],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
                set(gca, 'Xlim',[0 1000],'XTick',[0:200:1000])
set(gcf,'Position',[10 10 400 400])

figure(6)
plot(vmi0,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(vmi,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(f)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('V_m [m s^{-1}]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
        line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
            set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])



figure(7)
plot(zi0*1000.,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(zi*1000.,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(g)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('Z_{i,tot} [%%]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
%         line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%             set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])

figure(8)
plot(mui0,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(mui,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(h)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('mu [%%]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
%         line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%             set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])


figure(9)
plot(lami0,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(lami,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(i)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('lam [%%]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
%         line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%             set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])


figure(10)
plot(dm0*100.,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(dm*100.,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(i)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('Dm [cm]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
%         line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%             set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])


figure(11)
plot(ze0,z0*10.^-3,'Color','k','LineWidth',1)
hold on
plot(ze,z*10.^-3,'Color','k','Linestyle','--','LineWidth',3.5)
    ylabel('z [km]')
    title('(i)')
            set(gca, 'ylim',[0 1],'YTick',[0:0.5:1])
    xlabel('Ze [dBz]')
                    set(gca,'Fontsize',20,'FontName','Arial','FontWeight','bold')
%         line('XData',[0 8],'YData',[0.5 0.5],'Linewidth',2.5,'Color','k','LineStyle',':')
%             set(gca, 'Xlim',[0 8],'XTick',[0:2:8])
set(gcf,'Position',[10 10 400 400])