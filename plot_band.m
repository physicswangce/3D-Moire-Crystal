%band data generated by test_specturm
figure1 = figure;
tsize = 28;
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.35])
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
color = [253, 109,90;  254,180, 11;  109, 195, 84; 153, 68, 135; 81, 140,216; 68, 50, 149]/255;
% Create multiple lines using matrix input to plot
for p = 1 : 6
    plot(Specturm(:,p)/abs(t_),'LineWidth',1,'Parent',axes1,'color',[0 0 0]);
    
end

for k = 1 : 6
    plot([position(k+1) position(k+1)],[-18 18],'linestyle','--','LineWidth',0.2,'Parent',axes1,'color',[169 169 169]/255);
end
set(gca,'xcolor',[169 169 169]/255);
set(gca,'ycolor',[169 169 169]/255);
%plot(Specturm1,'LineWidth',1,'Marker','.','Parent',axes1);

%Set the remaining axes properties
set(gca,'TicklabelInterpreter','latex');
set(axes1,'FontSize',18,'GridAlpha',0.0, 'Linewidth',2 ,...
    'GridLineStyle','--','XGrid','on')
% set(gca,'XTick',position_1','XTickLabel',...
%     {'$\Gamma$','$$M$$','$$K$$','$$\Gamma$$','$$A$$','$$L$$','$$H$$','$$A$$'});
set(gca,'XTick',position_1','XTickLabel',{})
% set(gca,'YTick',[-10 0 10],'YTickLabel',...
%     {'$-10$','$0$','$10$'});
set(gca,'YTick',[-10 0 10],'YTickLabel',{})
set(gca,'Fontname','Times New Roman','Fontsize',0.1)
xlim([0 size(Kpath1,2)])
%ylabel('$$E/t$$','Fontname','Times New Roman','Fontsize',20)
ylabel('$E/t$','interpreter','latex','Fontsize',tsize,'position',[-38,0],'color',[0 0 0])
ylim([-18 18])

box on

set(gca,'position',[0.07,0.07,0.9,0.9]);
set(gca,'DataAspectRatio',[15 1 1])
po_=get(axes1,'Position');
indexlist = cell(8,1);
indexlist{1} = '$\Gamma$';indexlist{2} = '$$M$$';indexlist{3} = '$$K$$';indexlist{4} = '$$\Gamma$$';
indexlist{5} = '$$A$$';indexlist{6} = '$$L$$';indexlist{7} = '$$H$$';indexlist{8} = '$$A$$';

poslist = zeros(8,1);
poslist(1) = po_(1) - 0.015;
poslist(2) = po_(1) + 0.13;
for p = 3 : 8
    poslist(p) = poslist(1) + position_1(p)*(poslist(2)-poslist(1))/position_1(2);
end
poslist(2) = poslist(2)-0.01;
poslist(3)  = poslist(3)-0.002;
for p = 1 : 8
annotation(figure1,'textbox',...
    [poslist(p) po_(2)+0.035 0 0],...
    'String',{indexlist{p}},...
    'Interpreter','latex',...
    'FontSize',tsize,...
    'FitBoxToText','off',...
    'EdgeColor','none');

end

annotation(figure1,'textbox',...
    [0.006 0.35 0 0],...
    'String',{'-$10$'},...
    'Interpreter','latex',...
    'FontSize',tsize,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.038 0.59 0 0],...
    'String',{'$0$'},...
    'Interpreter','latex',...
    'FontSize',tsize,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.020 0.82 0 0],...
    'String',{'$10$'},...
    'Interpreter','latex',...
    'FontSize',tsize,...
    'FitBoxToText','off',...
    'EdgeColor','none');