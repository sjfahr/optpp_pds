%This parses the power further so that there are no repeats, redundant
%intervals.

function [Power_intervals]=power_parser_write_DF_array(power_log);

%Find where the power changes at all;
delta_P = (diff( power_log(:,5) )~=0) + ( diff(power_log(:,6) )~=0);  %Find which elements change from columns 5 and 6; then add the changes into one column
delta_P(1+1:end+1,:) = delta_P (1:end,:);  %add row back in because diff function eliminates the first one.
delta_P(1,:) = delta_P (2,:);
%At this point, delta_P lists all the times that columns 5 and 6 change

% Check to see if there is anything changing
if isempty ( find ( delta_P )) == 1;
    
    Power_intervals (1,1) = power_log ( end, 4 );
    Power_intervals (1,2) = power_log ( end, 6 );
    
    return
end

keep = find ( delta_P) ; % The +1 is to match the indexing (diff drops the length by 1)
on_off = zeros ( length( find ( delta_P ) )  , 1);

%This loop captures the on/off state of the power_log
for ii = 1 : length ( keep )
    
    on_off (ii) = power_log( keep(ii) , 5 );
    
end

clear ii

P = find (delta_P); %Column 1 of P records the times that the power changes

P(:,2) = power_log (P(:,1),6); %Use the times from column 1, P to record the corresponding powers


P(:,2) = P(:,2) * 15/100; %Convert % power to W power
P ( :,3 ) = on_off;

k_P (:,1)=P(:,1);
k_P (:,2)=P(:,2).*P(:,3);

times = k_P( : , 1 );
powers = k_P( :,2);

% Add a 0 power at the beginning.
powers = cat (1,0,powers);

times(end+1)=power_log(end,4);
%times(1) = []; A%%% prolly will cut
%times = round (times/2);

Power_intervals (:,1) = times;
Power_intervals (:,2) = powers;

Power_intervals_diff = diff (Power_intervals);

cut_indices = find ( not ( Power_intervals_diff ( :,2) ));

Power_intervals ( cut_indices , : )= [];

end