% plot_front.m
% 
% Matlab source codes
% 
% Plot the PF and PS of the test instances for CEC 2009 MOO Competition
%
% Usage: plot_front(), the result figures are in pf_figs folder
% 
% Please refer to the report for more information.
%

function plot_front()

plot_unconstraint();

plot_constraint();

end

%% ss
function plot_unconstraint()
PROBLEMS= ['UF1 '; 'UF2 '; 'UF3 '; 'UF4 '; 'UF5 '; 'UF6 '; 'UF7 '; 'UF8 '; 'UF9 '; 'UF10';];
DIMX    = [30 30 30 30 30 30 30 30 30 30];
NOP     = [100 400 100 100 100 100 100 5000 5000 5000]; 
PROPERTY= ['b-'; 'b-'; 'b-'; 'b-'; 'b.'; 'b.';'b-';'b.';'b.';'b.';];

for p=1:10
    [PF,PS] = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );
    RA      = xboundary( deblank(PROBLEMS(p,:)), DIMX(p) );  
    
    set(0,'units','centimeters');
    position=[0 0 17 6.5];
    h=figure;
    set(h,'PaperType','A4'); 
    set(h,'PaperUnits','centimeters'); 
    set(h,'paperpositionmode','auto');
    set(h,'PaperPosition',position);
    set(h,'units','centimeters');
    set(h,'position',position);
    hold off;

    subplot(1,2,1);
    if size(PF,1)== 3
        plot3(PF(1,:),PF(2,:),PF(3,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;
    else
        plot(PF(1,:),PF(2,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;
    end
        
    set(gca,'FontSize',8);
    xlabel('f1');ylabel('f2');
    title('Pareto front');
    xlim([0 1.2]); ylim([0 1.2]);  
    set(gca,'XTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});       
    set(gca,'YTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'YTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});             
    if size(PF,1)== 3
        set(gca,'ZTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'ZTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});  
        view([-45 10]);
    end
    grid off; box on;

    subplot(1,2,2);
    plot3(PS(1,:),PS(2,:),PS(3,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;  
    set(gca,'FontSize',8);
    xlabel('x1');ylabel('x2'); zlabel('x3');
    title('Pareto set');
    xlim(RA(1,:)); ylim(RA(2,:)); zlim(RA(3,:));            
    view([-45 10]); grid off; box on;

    f = sprintf('pf_figs/%s.eps',deblank(PROBLEMS(p,:)));
    saveas(h,f,'psc2');              
%      f = sprintf('pf_figs/%s.bmp',deblank(PROBLEMS(p,:)));
%      saveas(h,f);              

    close(h);
    set(0,'units','pixel');        
end
end

%% ss
function plot_constraint()
PROBLEMS= ['CF1 '; 'CF2 '; 'CF3 '; 'CF4 '; 'CF5 '; 'CF6 '; 'CF7 '; 'CF8 '; 'CF9 '; 'CF10';];
DIMX    = [10 10 10 10 10 10 10 10 10 10];
NOP     = [500 500 500 500 500 500 500 5000 5000 5000]; 
PROPERTY= ['b.'; 'b.'; 'b.'; 'b.'; 'b.'; 'b.';'b.';'b.';'b.';'b.';];

for p=1:10
    [PF,PS] = pareto( deblank(PROBLEMS(p,:)), NOP(p), DIMX(p) );
    
    set(0,'units','centimeters');
    position=[0 0 8 6.5];
    h=figure;
    set(h,'PaperType','A4'); 
    set(h,'PaperUnits','centimeters'); 
    set(h,'paperpositionmode','auto');
    set(h,'PaperPosition',position);
    set(h,'units','centimeters');
    set(h,'position',position);
    hold off;

    if size(PF,1)== 3
        plot3(PF(1,:),PF(2,:),PF(3,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;
    else
        plot(PF(1,:),PF(2,:),deblank(PROPERTY(p,:)),'MarkerSize',4); hold on;
    end
    set(gca,'FontSize',8);
    xlabel('f1');ylabel('f2');
    title('Pareto front');
    xlim([0 1.2]); ylim([0 1.2]);  
    set(gca,'XTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'XTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});       
    set(gca,'YTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
    set(gca,'YTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});             
    if size(PF,1)== 3
        set(gca,'ZTick',     [0 0.2 0.4 0.6 0.8 1.0 1.2]);
        set(gca,'ZTickLabel',{'0.0','0.2','0.4','0,6','0.8','1.0','1.2'});  
        view([-45 10]);
    end
    grid off; box on;

    f = sprintf('pf_figs/%s.eps',deblank(PROBLEMS(p,:)));
    saveas(h,f,'psc2');              
    close(h);
    set(0,'units','pixel');        
end
end
