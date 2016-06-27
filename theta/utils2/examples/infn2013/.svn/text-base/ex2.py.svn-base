
execfile("common.py")


# 2.a.

# The p-value for the hypothesis test with null hypothesis s=s0 versus s<s0
# for a counting experiment with known background b, in which the number of
# observed events n follows a Poisson distribution around s+b.
#def get_pvalue(s0, b, nobs):
    #TODO: implement
    
    

# for fixed values of b and nobs, scan through s0 from s0min to s0max using nscan points
# and print / plot the result.
def scan_s0_pvalue(b, nobs, s0min, s0max, nscan):
    delta_s = (s0max - s0min) / nscan
    pvals, svals = [], []
    for i in range(nscan): # i=0...nscan-1
        s0 = s0min + i*delta_s
        p = None #TODO: replace this line with proper p-value calculation and print the result!
        
        svals.append(s0)
        pvals.append(p)
    plot_xy(svals, pvals, 'p-vs-s.pdf')
        

# 2.b.

# return the maximum value for nmin such that the probability
# to observe n>=nmin for a Poisson
# with mean mu is >= pmin.
def find_nmin_poisson(pmin, mu):
    # note: this is a very inefficient implementation, but at least it should
    # be very transparent how it works (scipy.stats.poisson.ppf provides a more efficient
    # implementation).
    n = 1
    while poisson_p_ge(n, mu) > pmin: n+= 1
    n -= 1
    return n


# cl is the confidence level (=1-alpha)
def construct_belt(b, s0min, s0max, nscan, cl = 0.95):
    delta_s = (s0max - s0min) / nscan
    nmins, svals = [], []
    for i in range(nscan):
        s0 = s0min + i*delta_s
        nmin = None #TODO: replace this with the proper calculation
        
        svals.append(s0)
        nmins.append(nmin)
    plot_xy(svals, nmins, 'neyman_belt.pdf', ymin = 0, xlabel = 's', ylabel = 'n')
    



# 2.c.

# get the p-value for the hypothesis test with null hypothesis s=0
# using lower values of n as sign of deviation.
#def get_pvalue_bonly(b, nobs):
    #TODO: implement here!


def get_cls(s0, b, nobs):
    psb = get_pvalue(s0, b, nobs)
    pb = get_pvalue_bonly(b, nobs)
    return psb / pb
    
    
# TODO: implement scan_s0_clsvalue(b, nobs, s0min, s0max, nscan) here!
    



# 2.d.
# generate ntoy random numbers around mu0+-delta_mu by first throwing
# a random number for mu and then using this as Poisson mean.
# see (solutions to) exercise 1.c.
def generate_poisson_unc(mu0, delta_mu, ntoy):
    result = []
    for i in range(ntoy):
        theta = norm.rvs()
        mu = exp(math.log(1 + delta_mu/mu0) * theta) * mu0
        result.append(poisson.rvs(mu))
    return result

# count number of elements in l which are <=x (in analogy to "count_ge" in exercise 1):
def count_le(l, x):
    n=0
    for x0 in l:
        if x0 <= x: n+=1
    return n


#TODO: implement the routines below
#def get_pvalue_toys(s0, b0, delta_b,  nobs, ntoys):
#def get_pvalue_bonly_toys(b0, delta_b, nobs, ntoys):
#def get_cls_toys(s0, b0, delta_b, nobs, ntoys):
#def scan_s0_clsvalue_toys(b0, delta_b, nobs, s0min, s0max, nscan=10, ntoys = 1000):



# 2.e.

#model = build_shape_model()
#limit, limit_error = get_observed_cls_limit(model)
#print "limit for mu: %.3f +- %.3f" % (limit, limit_error)


