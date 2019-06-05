%% Spreading law analysis
%% Determines the spreading law for an oil droplet across a substrate
%% Elliot Goodwin

%% ----- List of variables -------------------------------
%        
%        t - time
%        r - radius of droplet at time t
%       dr - standard deviation of r
%        v - velocity of droplet at time t
%        h - height of droplet at time t
%        V - volume of droplet (constant, 7.6pL)
%    theta - angle between droplet and substrate at time t
%% -------------------------------------------------------

close all;
clear all;
format long;

%% IMPORT DATA

% import data set 1
data1 = importdata('drop2data1.txt');

% split the columns and append to vectors 't1' and 'r1'
t1 = data1(:,1);
r1 = data1(:,2);
clear data1;

% import data set 2
data2 = importdata('drop2data2.txt');
t2 = data2(:,1);
r2 = data2(:,2);
clear data2;

% find the mean and std of r
rtemp = [r1, r2];
ttemp = [t1, t2];

r = mean(rtemp,2);
dr = std(rtemp,0,2);
t = mean(ttemp,2);

clear rtemp, ttemp;

% find v
% velocity defined as dr/dt
% approximate gradient as deltar/deltat

for n = 1:size(r)-2;    
    v(n) = (r(n+2)-r(n))./(t(n+2)-t(n));
    dv(n) = ((dr(n+2).^2 + dr(n).^2)).^0.5 ./ (t(n+2) - t(n))';
end;
clear n;

% remove first and last values of t matrix (cannot define gradient
% at these points from given data)
t_v = t(2:end-1,:);

% define the function relating r, h, and V (assuming droplet is a spherical cap)
% solve function for h at each value of r
h = fsolve(@(h)(pi./6).*h.*(3.*r.^2 + h.^2) - 7600, r);
theta = (pi./2) - atan((r.^2 - h.^2)./(2.*h.*r));

% evaluate function to check for correct solutions
h_check = (pi./6).*h.*(3.*r.^2 + h.^2) - 7600;

% remove first and last values to agree with time data
theta_v = theta(2:end-1,:);

% compute uncertainty in theta
a1 = (dr./(h.^2+r.^2));
a2 = (7600./(pi.*h.^2)+(h./3)).^-2;
dtheta = a1.*(2.*h.^2 + ((r.^2)./2).*a2).^0.5;



%%  FITTING
% quadratic fit; non-uniformly distributed data
[p2, s2] = polyfit(theta_v, v', 2);
y2 = p2(1).*(theta_v).^2 + p2(2).*(theta_v) + p2(3);

% find estimates of error at each point
err_v = ones(size(dv)).*(mean(dv));
% chi square
temp2 = ((v' - y2)./dv');
chisq2 = sum(temp2.^2);
chisq21 = (((s2.normr).^2)./(err_v(1).^2));

% find residuals
error_v2 = (s2.normr)./sqrt(s2.df);
sigma_v2 = error_v2*ones(size(theta_v));
vfit2 = polyval(p2, theta_v, 2);

% find errors on fit parameters
% covariance matrix
covm2 = error_v2*inv(s2.R)*inv(s2.R)';
% extract errors on parameters
err_p21 = sqrt(covm2(1,1));
err_p22 = sqrt(covm2(2,2));
err_p23 = sqrt(covm2(3,3));

% plot quadratic fit
figure(1);
errorbar(theta_v, v, dv, 'bo'); hold on;
plot(theta_v, y2, 'r-');  hold on;
% add axes labels
xlabel('\theta [rad]'); hold on;
ylabel(' velocity' );  hold off;
% plot residuals
figure(2)
errorbar(theta_v, (vfit2 - v'), sigma_v2, 'bo');    hold on;
xlabel('Fitted value'); hold on;
ylabel('Residual');  hold off;

% cubic fit; non-uniformly distributed data
[p3, s3] = polyfit(theta_v, v', 3);
y3 = p3(1).*(theta_v).^3 + p3(2).*(theta_v).^2 + p3(3).*(theta_v) + p3(4);

% chi square
temp3 = ((v' - y3)./dv');
chisq3 = sum(temp3.^2);
chisq31 = (((s3.normr).^2)./(err_v(1).^2));

figure(3);
errorbar(theta_v, v, dv, 'ro'); hold on;
plot(theta_v, y3);  hold on;
xlabel('\theta [rad]'); hold on;
ylabel('velocity [\micro ms^-1]');  hold off;

% find residuals
error_v3 = (s3.normr)./sqrt(s3.df);
sigma_v3 = error_v3*ones(size(theta_v));
vfit3 = polyval(p3, theta_v, 3);
% plot residuals
figure(4)
errorbar(theta_v, (vfit3 - v'), sigma_v3, 'bo');    hold on;
xlabel('Fitted value'); hold on;
ylabel('Residual');  hold off;

% find errors on fit parameters
% covariance matrix
covm3 = error_v3*inv(s3.R)*inv(s3.R)';
% extract errors on parameters
err_p31 = sqrt(covm3(1,1));
err_p32 = sqrt(covm3(2,2));
err_p33 = sqrt(covm3(3,3));
err_p34 = sqrt(covm3(4,4));



%% FITTING (Uniform data)
% spline data
xx = [theta_v(end):0.001:theta_v(1)];
yy = spline(theta_v, v, xx);


% quadratic fit; uniform data
[p4, s4] = polyfit(xx, yy, 2);
y4 = p4(1).*(xx).^2 + p4(2).*(xx) + p4(3);

% find estimates of error at each point on splined data
err_vsplined = ones(size(xx)).*(mean(dv));
% chi square
temp4 = ((yy - y4)./err_vsplined);
chisq4 = sum(temp4.^2);
% compare to estimate from polyfit data
chisq41 = ((s4.normr).^2)./(err_vsplined(1).^2);

%divide by s.df (degrees of freedom) to find reduced chi square
red_chisq4 = (sum(temp4.^2))./(s4.df);
red_chisq41 = (((s4.normr).^2)./(err_vsplined(1).^2))./(s4.df);

% show plot
figure(5);
errorbar(theta_v, v, err_v, 'ro'); hold on;
plot(xx, y4);  hold on;
xlabel('\theta [rad]'); hold on;
ylabel('velocity [micrometers/s]');  hold off;

% find residuals
error_v4 = (s4.normr)./sqrt(s4.df);
sigma_v4 = error_v4*ones(size(theta_v));
vfit4 = polyval(p4, theta_v, 2);
% plot residuals
figure(6);
errorbar(theta_v, (vfit4 - v'), sigma_v4, 'bo');    hold on;
xlabel('Fitted value'); hold on;
ylabel('Residual');  hold off;

% find errors on fit parameters
% covariance matrix
covm4 = error_v4*inv(s4.R)*inv(s4.R)';
% extract errors on parameters
err_p41 = sqrt(covm4(1,1));
err_p42 = sqrt(covm4(2,2));
err_p43 = sqrt(covm4(3,3));



% cubic fit; uniform data
[p5, s5] = polyfit(xx, yy, 3);
y5 = p5(1).*(xx).^3 + p5(2).*(xx).^2 + p5(3).*(xx) + p5(4);

% chi square
temp5 = ((yy - y5)./(err_vsplined));
chisq5 = sum(temp5.^2);
% compare to estimate from polyfit data
chisq51 = ((s5.normr).^2)./(err_vsplined(1).^2);

%divide by s.df (degrees of freedom) to find reduced chi square
red_chisq5 = (sum(temp5.^2))./(s5.df);
red_chisq51 = (((s5.normr).^2)./(err_vsplined(1).^2))./(s5.df);

% show plot
figure(7);
errorbar(theta_v, v, err_v, 'ro'); hold on;
plot(xx, y5);  hold on;
xlabel('\theta [rad]'); hold on;
ylabel('velocity');  hold off;

% find residuals
error_v5 = (s5.normr)./sqrt(s5.df);
sigma_v5 = error_v5*ones(size(theta_v));
vfit5 = polyval(p5, theta_v, 3);
% plot residuals
figure(8)
errorbar(theta_v, (vfit5 - v'), sigma_v5, 'bo');    hold on;
xlabel('Fitted value'); hold on;
ylabel('Residual');  hold off;

% find errors on fit parameters
% covariance matrix
covm5 = error_v5*inv(s5.R)*inv(s5.R)';
% extract errors on parameters
err_p51 = sqrt(covm5(1,1));
err_p52 = sqrt(covm5(2,2));
err_p53 = sqrt(covm5(3,3));
err_p54 = sqrt(covm5(4,4));

% save figures in png format
saveas(figure(1),'quad_fit','png');
saveas(figure(2),'quad_fit_residuals','png');
saveas(figure(3),'cube_fit','png');
saveas(figure(4),'cube_fit_residuals','png');
saveas(figure(5),'quad_spline','png');
saveas(figure(6),'quad_spline_residuals','png');
saveas(figure(7),'cube_spline','png');
saveas(figure(8),'cube_spline_residuals','png');

close all;
