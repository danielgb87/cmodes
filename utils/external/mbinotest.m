
function p = mbinotest(nheads,ntosses,H0p,tails,direction)
% function p = mbinotest(nheads,ntosses,H0p,tails,direction)
% conducts a binomial one- or two-tailed test
% assessing the probability of obtaining the observed pattern of a Bernoulli experiment under the null hypothesis
% for simplicity, we will use the classical experiment of tossing a biased coin for illustration
%
% INPUT ARGUMENTS
% nheads          number of heads obtained after tossing the coin a specified number of times
% ntosses         how many coin tosses were performed
% H0p             what fraction of tosses did we expect to yield heads under the null hypothesis (0.5 for a fair coin)
% tails           do we want a one- or a two-sided p-value
% direction       optional argument, only used when tails==1
%                 'more' or 'less' allocates the rejection region to the upper and lower range of outcomes, respectively
%
% OUTPUT ARGUMENT
% p     probability to observe nheads or more/less with ntosses under the null hypothesis of H0p
%
% EXAMPLES
% Suppose we have flipped a coin 80 times and obtained heads 45 times. We think the coin is far (H0p = 0.5).
% p = mbinotest(45,80,0.5,2);
% We obtain p = 0.31431, an indication that the coin is fair.
%
% Now suppose we think the coin is biased towards landing heads (i.e. H0p>0.5).
% p = binotest(45,80,0.5,1,'more');
% p = 0.15715, so no strong support for our suspicion. Note that for a 2-tailed test, the p-value had been twice is high. Also note that,
% had we expected less heads with same results: p = binotest(45,80,0.5,1,'less'), our p-value would have equaled 0.8907 instead. This is
% because the p-value is an integral over all outcomes that are "as extreme or more extreme" than the outcome actually observed.
%
% August 2016 by Maik C. Stuettgen, University Medical Center Mainz, Germany
%% input check
if rem(nheads,1) ~=0 || rem(ntosses,1)~=0 || nheads<1 || ntosses<1
  error('arguments ''nheads'' and ''ntosses'' must both be integers')
elseif nheads>ntosses
  error('argument ''nheads'' must be smaller than ''ntosses''')
elseif H0p<=0 || H0p>=1
  error('argument ''H0p'' must be a value between 0 and 1')
elseif ~ismember(tails,[1 2])
  error('argument ''tails'' must be either 1 or 2')
elseif exist('direction','var') && all([~strcmp(direction,'more') ~strcmp(direction,'less')])
  error('argument ''direction'' must read either ''more'' or ''less''')
elseif tails==1 && ~exist('direction','var')
  error('argument ''direction'' not specified')
end
%% the works
if tails==1   % this could be shortened but I like it because its nicely ordered
  if nheads>ntosses*H0p && strcmp(direction,'more')         % nheads in predicted direction 'more'
    p = 1-binocdf(nheads-1,ntosses,H0p);
  elseif (nheads<ntosses*H0p && strcmp(direction,'less'))   % nheads in predicted direction 'less'
    p = binocdf(nheads,ntosses,H0p);
  elseif (nheads>ntosses*H0p && strcmp(direction,'less'))   % nheads not in predicted direction 'less'
    p = binocdf(nheads,ntosses,H0p);
  elseif nheads<ntosses*H0p && strcmp(direction,'more')     % nheads not in predicted direction 'more'
    p = 1-binocdf(nheads-1,ntosses,H0p);
  elseif nheads==ntosses*H0p                                % nheads exactly equals expected value under the null hypothesis
    p = binocdf(nheads,ntosses,H0p);
  end
elseif tails==2
  if nheads>ntosses*H0p
    p = 2*(1-binocdf(nheads-1,ntosses,H0p));
  elseif nheads<=ntosses*H0p
    p = 2*binocdf(nheads,ntosses,H0p);
  end
end
end
