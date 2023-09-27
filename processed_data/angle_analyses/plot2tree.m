function y = plot2tree(f_name)

fid = fopen(f_name);
C = textscan(fid, '%f %f %f %f');
fclose(fid);

L= numel(C{1});

figure
hold on

for ii = 1:50:L
    plot([C{1}(ii) C{3}(ii)],[C{2}(ii) C{4}(ii)],'-o',...
     'LineWidth',3, 'MarkerSize',4,'MarkerEdgeColor','k','Color','blue')
    plot([C{1}(ii+1) C{3}(ii+1)],[C{2}(ii+1) C{4}(ii+1)],'-o',...
     'LineWidth',3, 'MarkerSize',4,'MarkerEdgeColor','k','Color','red')
end

y=0;