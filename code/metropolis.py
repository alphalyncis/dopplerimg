import numpy as np
from scipy.special import comb, factorial
from numpy import log

n_det = 10


# All the probability functions here are actually log probabilities.

# Priors
def lam_prior(lam):
    return 0.0
    
def f_prior(f):
    return -0.5*(f - 0.5)**2/0.05**2



# Likelihood components
def log_factorial(n):
    return sum(log(i) for i in range(1,n+1))

def n_emit_prob(n_emit, lam):
    return log(lam)*n_emit - lam -  log_factorial(n_emit)

def n_det_prob(n_emit, f):
    return log(comb(n_emit, n_det)) +  n_det*log(f) +  (n_emit-n_det)*log(1-f)

# Total posterior = prior + likelihood
def posterior(n_emit, f, lam):
    prior = lam_prior(lam) + f_prior(f)
    likelihood = n_emit_prob(n_emit, lam) + n_det_prob(n_emit,f)
    return prior + likelihood


# Sampling
def proposal(n_emit, f, lam):
    n_emit = n_emit + np.random.choice([-4,-3,-2,-1,0,1,2,3,4])
    f = f + np.random.normal()*0.05
    lam = lam + np.random.normal()*2.
    return n_emit, f, lam

def metropolis(N):
    #Fixed starting point appropriate for this problem
    x = [20, 0.5, 20.0]
    logp = posterior(*x)
    chain = []
    n = 0
    # Main loop of the algorithm
    for i in range(N):
        # The * thing converts a list like x to a sequence of arguments.
        # i.e. if x = [1,2,3] then f(*x) = f(1,2,3)
        x_new = proposal(*x)
        logp_new = posterior(*x_new)
        # Metropolis-Hastings acceptance test
        if (logp_new-logp) > log(np.random.uniform()):
            x = x_new
            logp = logp_new
            n += 1
        # Store result
        chain.append(x)
        #Print out current acceptance rate
        if i>0 and (i%1000 == 0):
            print(n*100./i)
    return chain

# Save results to disc
def save(chain, filename):
    f = open(filename,'w')
    for row in chain:
        f.write("{}  {}  {}\n".format(*row))
    f.close()


if __name__ == '__main__':
    chain = metropolis(100000)
    save(chain, 'chain.txt')