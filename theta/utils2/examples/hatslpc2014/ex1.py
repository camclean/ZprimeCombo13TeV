execfile("common.py")


# 1.a.
def ex1a():
    #b = 5.2
    b = 0.4
    print "p-values for a counting experiment with b=", b
    for nobs in [1,2,3]: #[5,6,7,8,10,15,20]:
        p = poisson_p_ge(nobs, b)
        Z = p_to_Z(p)
        print "for nobs=%d: p=%.3g; Z=%.2f" % (nobs, p, Z)
        #Z_a = (nobs - b)/sqrt(b)
        #print "for nobs=%d: p=%.3g; Z=%.2f; Z_a=%.2f" % (nobs, p, Z, Z_a)

# calculate the approximate Z-value based on Wilks' Theorem for a counting experiment
# with known b (for 1.a. iii.)
def calc_Z_w(n, b):
    s_hat = max([0.0, n - b])
    return sqrt(2 * ( n * log((s_hat + b) / b) - s_hat))

ex1a()





# 1.b.

def generate_poisson(mu, ntoy):
    return poisson.rvs(mu, size = ntoy)

# I. print some generated Poisson numbers:
#data = generate_poisson(5.0, 10)
#print data


# II. plot the distribution of generated Poisson data; "poisson5.pdf" will be created
#    in the current working directory you call theta-auto from.
#data = generate_poisson(5.0, 1000)
#plot_histogram(data, 'poisson5.pdf')


# return the number of elements in the list l which are >= x
def count_ge(l, x):
    n=0
    for x0 in l:
        if x0 >= x: n+=1
    return n
    # equivalent implementation: return sum([1 for x0 in l if x0 >= x])


def get_pvalue(b, nobs, ntoy = 1000):
    # generate an ensemble of values of n for background only:
    bonly_ns = generate_poisson(b, ntoy)
    # count the number of toys for which n >= nobs:
    ntoy_ge_nobs = count_ge(bonly_ns, nobs)
    # the estimated p-value is the fraction of toys with n >= nobs:
    return ntoy_ge_nobs * 1.0 / ntoy

# III. get the p-value for b=5.2 and nobs = 8:
#p = get_pvalue(5.2, 8)
#print "p=%.3f; Z=%.3f" % (p, p_to_Z(p))




# 1.c.

def generate_poisson_unc(mu0, delta_mu, ntoy):
    result = []
    for i in range(ntoy):
        # generate Poisson mean mu according to a normal distribution with mean mu0 and width delta_mu
        mu = norm.rvs(mu0, delta_mu)
        # generate a Poisson random number around mean mu and append it to the result list:
        result.append(poisson.rvs(mu))
    return result

# Using get_pvalue as starting point, write a method called "get_pvalue_syst" for calculating the p-value with a
# background uncertainty.


# 1.d.
def ex1d():
    # In theta, the statistical model is an instance of the class "Model",
    # which contains all information for other routines -- including the observed data.
    # This is the shape model introduced in the lecture:
    model = build_shape_model(signal_mass = 500.)
    # get the test statistic values for background-only:
    ntoy = 5000
    t_bkg = get_bkg_t(model, ntoy)
    tobs = get_data_t(model)
    # note: you can use plot_histogram as in 1.b. to visualize the test statistic distribution
    # count the number of toys for which t >= tobs:
    n = count_ge(t_bkg, tobs)
    p = n * 1.0 / ntoy
    Z = p_to_Z(p)
    print "p = %.3g; Z = %.3g" % (p, Z)
    
    
#ex1d()

