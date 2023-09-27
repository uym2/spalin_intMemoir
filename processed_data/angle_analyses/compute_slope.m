function slope = compute_slope(f_name)

fid = fopen(f_name);
C = textscan(fid, '%f %f %f %f');
fclose(fid);

x_p = C{1}(1:2:end);
y_p = C{2}(1:2:end);
x_c1 = C{3}(1:2:end);
y_c1 = C{4}(1:2:end);
x_c2 = C{3}(2:2:end);
y_c2 = C{4}(2:2:end);

v = [x_c1-x_c2,y_c1-y_c2];
mean(abs(v(v(:,1)==0,2)))
mean(abs(v(:,2)))
mean(abs(v(v(:,2)==0,1)))
mean(abs(v(:,1)))
slope = round(atan(v(:,2)./v(:,1))/pi*180);