clear all;

N=100;
GoldenRatioS=1.5+0.5*sqrt(5);   % the square of the golden ratio
alphaGA=2*pi/GoldenRatioS;      % the golden angle in radians
alpha=alphaGA;
%  alpha = 45/360*pi;
%
VogelArrayXY=zeros(2,N);
VogelArrayRT=VogelArrayXY;
min_dist = 400;
a=1/1.6*min_dist;
for ip=1:N;
    rp=sqrt(ip)*a;
    alphap=alpha*ip;
    VogelArrayRT(1,ip)=rp;
    VogelArrayRT(2,ip)=alphap;
    VogelArrayXY(1,ip)=rp*cos(alphap);
    VogelArrayXY(2,ip)=rp*sin(alphap);
end

% figure; plot(VogelArrayXY(1,:),VogelArrayXY(2,:),'.'); axis equal;
% min_dist = 100000;
% for ip=1:N
%     for secip = 1:N
%         
%         if (ip==secip)
%             dista(ip,secip) = 10000; %%#
%         else
%             dista(ip,secip) = sqrt((VogelArrayXY(1,ip)-VogelArrayXY(1,secip)).^2+ (VogelArrayXY(2,ip)-VogelArrayXY(2,secip)).^2);
%         end
%         if (dista(ip,secip)<min_dist)
%             
%             min_dist = dista(ip,secip);
%             a_in = ip;
%             b_in = secip;
%         end
%     end
%     
%     min_dist_val_ind(1,ip) = min(dista(ip,:));
%     min_dist_val_ind(2,ip) = find(min_dist_val_ind(1,ip)==dista(ip,:));
%     min_dist_val_ind(3,ip) = ip;
% end
%  
% min_dist_val_ind_sorted = sort(min_dist_val_ind,2);
% min_dist_val_ind_sorted = min_dist_val_ind_sorted';
%    min_dist
%    a_in  
%    b_in
%    dista_sort = sort(dista(:));
%   dista_sort(dista_sort==0) = [];
% 
% 
% figure; plot(VogelArrayXY(1,:),VogelArrayXY(2,:),'.'); axis equal; 
% hold on; 
% plot(VogelArrayXY(1,a_in),VogelArrayXY(2,a_in),'or');
% plot(VogelArrayXY(1,b_in),VogelArrayXY(2,b_in),'og');