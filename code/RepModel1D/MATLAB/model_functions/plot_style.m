function plot_style(Xlabel,Ylabel,Title,Wsize,Wloc)
if ~exist('Wsize','var')
    Wsize = [1000 500]; 
end
if ~exist('Wloc','var')
    Wloc = [50 50];
end

set(gcf,'Renderer', 'painters', 'Position', [Wloc Wsize])
if ~isempty(Xlabel);    xlabel(Xlabel,"Interpreter","latex"); end
if ~isempty(Ylabel);    ylabel(Ylabel,"Interpreter","latex"); end
if ~isempty(Title);     title(Title,"Interpreter","latex"); end
set(gca,'fontname','times','Layer','top'); box on; grid on;
set(gcf,'color','w'); 
end
