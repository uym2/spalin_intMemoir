function y = plot3tree(f_name)

fid = fopen(f_name);
C = textscan(fid, '%f %f %f %f %f %f');
fclose(fid);

L= numel(C{1});

figure
hold on

x_all = [C{1}; C{4}];
x_range = max(x_all) -min(x_all);
x_lim = [min(x_all)-0.1*x_range, max(x_all)+0.1*x_range];

y_all = [C{2}; C{5}];
y_range = max(y_all) -min(y_all);
y_lim = [min(y_all)-0.1*y_range, max(y_all)+0.1*y_range];

[X, Y] = meshgrid(x_lim,y_lim);

Z = -216*ones(size(X));
surf(X,Y,Z,'FaceColor','c')

x1_lim = [min(x_all)+3*x_range/16, max(x_all)-3*x_range/16];
y1_lim = [min(y_all)+3*y_range/16, max(y_all)-3*y_range/16];
[X1,Y1] = meshgrid(x1_lim,y1_lim);
Z1 = -160*ones(size(X));
surf(X1,Y1,Z1,'FaceColor','c')

for ii = 1:L
    plot3([C{1}(ii) C{4}(ii)],[C{2}(ii) C{5}(ii)],[-C{3}(ii) -C{6}(ii)],'-o',...
        'LineWidth',3, 'MarkerSize',4,'MarkerEdgeColor','k')
end

% set(gca, 'Xtick',[])
% set(gca, 'Ytick',[])
set(gca, 'Ztick',[])
set(gca, 'CameraPosition', [0 0 0]);


y=0;