function theta = compute_theta(f_name)

fid = fopen(f_name);
C = textscan(fid, '%f %f %f %f');
fclose(fid);

x_p = C{1}(1:2:end);
y_p = C{2}(1:2:end);
x_c1 = C{3}(1:2:end);
y_c1 = C{4}(1:2:end);
x_c2 = C{3}(2:2:end);
y_c2 = C{4}(2:2:end);

v1 = [x_c1-x_p,y_c1-y_p];
v2 = [x_c2-x_p,y_c2-y_p];

cos_theta = sum(v1.*v2,2)./(sqrt(sum(v1.*v1,2)).*sqrt(sum(v2.*v2,2)));
%cos_theta = cos_theta(~isnan(cos_theta));
theta = acos(cos_theta)/pi*180;