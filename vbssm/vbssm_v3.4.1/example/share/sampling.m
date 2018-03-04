function sampled_data = sampling(data,t,t_wds,samples)
tstart = t_wds(1);
tend = t_wds(2);
startind = min( find( t>=tstart ) );
endind = min( find( t>=tend ) );

% number of timepoints
numtp = endind - startind;
interval = floor(numtp/(samples));
sampled_data = data(startind:interval:endind,11:20,:);
