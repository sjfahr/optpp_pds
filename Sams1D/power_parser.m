%This parses the power further

function [P,unique_P,delta_P]=power_parser(power_log);

%Find where the power changes at all;
delta_P=(diff(power_log(:,5))~=0)+(diff(power_log(:,6))~=0);  %Find which elements change from columns 5 and 6; then add the changes into one column
delta_P(1+1:end+1,:) = delta_P(1:end,:);  %add row back in because diff function eliminates the first one.
delta_P(1,:)=delta_P(2,:); 
%At this point, delta_P lists all the times that columns 5 and 6 change

delta_P=power_log(:,5).*delta_P;  %Only keeps power changes while laser is on (power setting can change while the laser is off, and I don't care about those changes)

P=find(delta_P); %Column 1 of P records the times that the power changes

P(:,2)=power_log(P(:,1),6); %Use the times from column 1, P to record the corresponding powers


P(:,2)=P(:,2)*15/100; %Convert % power to W power
unique_P=unique(P(:,2));
%At this point, P column 1 is time, column 2 is the power in W
% j=1;
% 
% if unique_P(i)==P(i,2)
%     unique_P(i,j)=j
%     j+1;

end