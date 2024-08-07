% The following is Lee and Valiant's toy implementation of their estimator, provided
% with their permission.

%% k experiment trials of n samples each
logdelta = 3;
% failure probability is e^{-logdelta}, the following code assumes logdelta
% divides n, so that we don't run into issues with the simple
% median-of-means implementation.
n = 99;
k = 1000000;

%% Example samples from the Student t distribution
nu = 3;
samples = random('T', nu, k, n);

%% MoM
groups = reshape(samples, k, logdelta, []);
group_means = mean(groups,3);
kappa = median(group_means,2);

%% actual estimator
dev = (samples-kappa);
sqdist_sorted = sort(dev.^2,2);
alphas_cand = (logdelta/3 - repmat(((n-1):-1:0),k,1))./cumsum(sqdist_sorted,2);
mask = ((alphas_cand.*[sqdist_sorted(:,2:end) inf(k, 1)]) >= 1) & (alphas_cand.*sqdist_sorted <= 1); % There might be slight numerical issues with this, but experiments have been ok
mask(:,n) = mask(:,n)+1*(sum(mask,2) < 1); % take the last alpha if too few samples to throw out
alpha = alphas_cand(mask == 1);
est = kappa + 1/n*sum(dev.*(1-min(alpha.*dev.^2,1)),2);

% You can in fact re-run the above section multiple times, reinitializing
% kappa as the previous estimate, and run the estimator until convergence.
% Our paper shows that even in 1 step, we get 1+o(1)-tightness, though
% it's possible that the fixpoint has a better o(1) term.