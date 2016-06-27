execfile("common.py")
execfile("ex1.py")

# for 1.a. and 1.b., the solutions are in ex1.py already, just comment in the
# right parts and change some numbers!


# 1.c.: fixed generate_poisson_unc:
def generate_poisson_unc(mu0, delta_mu, ntoy):
    result = []
    for i in range(ntoy):
        mu = norm.rvs(mu0, delta_mu)
        # The problem is here: mu can be <= 0 and a Poisson with mean
        # <= 0 does not make sense. A simple "solution" is to use a truncated normal distribution,
        # i.e. re-generate until we have a value mu >0:
        #while mu <=0: mu = norm.rvs(mu0, delta_mu)
        
        # Another (and better) possible fix is not to use a normal
        # prior on b but a log-normal prior. As discussed in the lecture, this
        # can be implemented by generating a normally distributed number
        # theta with mean 0 and standard deviuation 1 and using it in an 
        # exponential scale factor for mu0:
        
        theta = norm.rvs()
        mu = exp(math.log(1 + delta_mu/mu0) * theta) * mu0
        
        # generate a Poisson random number with mean mu:
        result.append(poisson.rvs(mu))
    return result
    
    
# 1.c.: implementation of get_pvalue_syst
def get_pvalue_syst(b, bunc, nobs, ntoy = 1000):
    # generate an ensemble of values of n for background only:
    bonly_ns = generate_poisson_unc(b, bunc, ntoy)
    # count the number of toys for which n >= nobs:
    ntoy_ge_nobs = count_ge(bonly_ns, nobs)
    # the estimated p-value is the fraction of toys with n >= nobs:
    return ntoy_ge_nobs * 1.0 / ntoy
    

# 1.c. code to answer the question on p-values, Z-values and comparison with
# noral approximation:
def ex1c_question():
    b = 10000.
    bunc = 50.
    nobs = 10200
    ntoy = 5000
    p = get_pvalue_syst(b, bunc, nobs, ntoy)
    Z = p_to_Z(p)
    Z_a = (nobs - b) / sqrt(b + bunc**2)
    
    # implementation of uncertainty on p: use normal approximation
    # for binomial error on p:
    delta_p = sqrt(p*(1-p)/ntoy)
    Zup = p_to_Z(p+delta_p)
    Zdown = p_to_Z(p-delta_p)

    print "Z_toy=%.2f (%.2f--%.2f); Z_a=%.2f" % (Z, Zup, Zdown, Z_a)
    

#ex1c_question()

# 1.d.

def ex1d_question():
    ntoy = 5000
    # 1. no uncertainty (see ex1.py):
    model = build_shape_model(signal_mass = 500.)
    t_bkg = get_bkg_t(model, ntoy)
    tobs = get_data_t(model)
    n = count_ge(t_bkg, tobs)
    p = n * 1.0 / ntoy
    Z = p_to_Z(p)
    print "no uncertainties:"
    print "p = %.3g; Z = %.3g" % (p, Z)
    
    # 2. with both uncertainties:
    model = build_shape_model(signal_mass = 500., b_rate_unc = 0.1, b_shape_unc = True)
    t_bkg = get_bkg_t(model, ntoy)
    tobs = get_data_t(model)
    n = count_ge(t_bkg, tobs)
    p = n * 1.0 / ntoy
    Z = p_to_Z(p)
    print "both rate and shape uncertainties:"
    print "p = %.3g; Z = %.3g" % (p, Z)
    
    
#ex1d_question()



