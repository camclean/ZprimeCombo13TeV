
execfile("common.py")
execfile("ex2.py")

# 2.a.

# The p-value for the hypothesis test with null hypothesis s=s0 versus s<s0
# for a counting experiment with known background b, in which the number of
# observed events n follows a Poisson distribution around s+b.
def get_pvalue(s0, b, nobs):
    # here, the p-value is the probability to observe <= n events
    # in a Poisson distribution with mean s0 + b:
    return poisson_p_le(nobs, s0+b)

# code to answer the question:
def ex2a_question_i():
    print "p-values for nobs=6, b=5.2 for s=s0 versus s<s0:"
    for s0 in [0, 1, 5, 10]:
        p = get_pvalue(s0, 5.2, 6)
        print "p for s0=%.2f: %.3g" % (s0, p)
    #
    # output:
    #p-values for nobs=6, b=5.2 for s=s0 versus s<s0:
    #p for s0=0.00: 0.732
    #p for s0=1.00: 0.574
    #p for s0=5.00: 0.118
    #p for s0=10.00: 0.00672
    #
    #Answer to ii: s0=0,1,5 are part of the interval as p>alpha=0.05, while s0=10 can be excluded
    #  and is therefore not part of the confidence interval for s.

#ex2a_question_i()


def scan_s0_pvalue(b, nobs, s0min, s0max, nscan):
    delta_s = (s0max - s0min) / nscan
    pvals, svals = [], []
    for i in range(nscan): # i=0...nscan-1
        s0 = s0min + i*delta_s
        p = get_pvalue(s0, b, nobs)
        print "p for s0=%.2f: %.3g" % (s0, p)
        svals.append(s0)
        pvals.append(p)
    plot_xy(svals, pvals, 'p-vs-s.pdf')

def ex2a_question_iii():
    # from ii, we know that the interval is somewhere between 5.0 and 10.0, so
    # scan this range more accurately:
    scan_s0_pvalue(5.2, 6, 5.0, 10.0, 100)
    
    # output:
    # ...
    #p for s0=6.55: 0.0526
    #p for s0=6.60: 0.0512
    #p for s0=6.65: 0.0498
    #p for s0=6.70: 0.0484
    #p for s0=6.75: 0.0471
    #p for s0=6.80: 0.0458
    #...
    #
    # -> 95% C.L. upper limit for s is about 6.65.

#ex2a_question_iii()

def ex2a_question_iv():
    scan_s0_pvalue(0, 0, 0.0, 5.0, 100)
    # -> answer: The 95% C.L. upper limit is s=3.0.

#ex2a_question_iv()

#2.b.
        
# cl is the confidence level (=1-alpha)
def construct_belt(b, s0min, s0max, nscan, cl = 0.95):
    delta_s = (s0max - s0min) / nscan
    nmins, svals = [], []
    for i in range(nscan):
        s0 = s0min + i*delta_s
        # for given s0, the poisson mean is mu = b+s0, and we want
        # to choose nmin such that the probability for n>=nmin is the confidence level cl:
        nmin = find_nmin_poisson(cl, s0 + b)
        svals.append(s0)
        nmins.append(nmin)
    plot_xy(svals, nmins, 'neyman_belt.pdf', ymin = 0, xlabel = 's', ylabel = 'n')
    
    # answer to question:
    # - 95% C.L. limit for s for n=6 is about 6.5 (see 2.a. for more accurate value)
    # - 95% C.L. limit for s for n=10 is about 12
    # - 95% C.L. limit for s for n=1 is 0.0, i.e. an empty interval
        

#construct_belt(5.2, 0.0, 15.0, 1000)


# 2.c.

def get_pvalue_bonly(b, nobs):
    return poisson_p_le(nobs, b)


def scan_s0_clsvalue(b, nobs, s0min, s0max, nscan):
    delta_s = (s0max - s0min) / nscan
    clsvals, svals = [], []
    for i in range(nscan): # i=0...nscan-1
        s0 = s0min + i*delta_s
        cls = get_cls(s0, b, nobs)
        # note: we only print value close to the threshold alpha=0.05,
        # so we are not overwhelmed by output on the console:
        if abs(cls-0.05) < 0.001: 
            print "CLs for s0=%.2f: %.3g" % (s0, cls)
        svals.append(s0)
        clsvals.append(cls)
    plot_xy(svals, clsvals, 'cls-vs-s-%d.pdf' % nobs)
    
    
def ex2c_question():
    for nobs in [6, 10, 1]:
        print "** nobs = %d" % nobs
        scan_s0_clsvalue(5.2, nobs, 0.0, 20.0, 1000)
    # output:
    #** nobs = 6
    #CLs for s0=7.18: 0.0505
    #CLs for s0=7.20: 0.05
    #CLs for s0=7.22: 0.0494
    #** nobs = 10
    #CLs for s0=11.76: 0.051
    #CLs for s0=11.78: 0.0505
    #CLs for s0=11.80: 0.05
    #CLs for s0=11.82: 0.0495
    #CLs for s0=11.84: 0.0491
    #** nobs = 1
    #CLs for s0=3.42: 0.0508
    #CLs for s0=3.44: 0.0499
    #
    #answer:
    # the limits are 7.2 for nobs=6, 11.8 for nobs=10 and 3.4 for nobs=1
   
#ex2c_question()


# 2.d.

def get_pvalue_toys(s0, b0, delta_b,  nobs, ntoys):
    ns = generate_poisson_unc(s0+b0, delta_b, ntoys)
    p = count_le(ns, nobs) * 1.0 / ntoys
    return p
    
def get_pvalue_bonly_toys(b0, delta_b, nobs, ntoys):
    ns = generate_poisson_unc(b0, delta_b, ntoys)
    p = count_le(ns, nobs) * 1.0 / ntoys
    return p
    
# note: here we also calculate the uncertainty on the CLs value
#  from the limited amount of toys. The return value is the tuple
# cls, delta_cls
def get_cls_toys(s0, b0, delta_b, nobs, ntoys):
    psb = get_pvalue_toys(s0, b0, delta_b, nobs, ntoys)
    pb = get_pvalue_bonly_toys(b0, delta_b, nobs, ntoys)
    delta_psb = sqrt(psb*(1-psb)/ntoys)
    delta_pb = sqrt(pb*(1-pb)/ntoys)
    cls = psb / pb
    # by the "usual" Gaussian error propagation,
    # the relative error on cls = psb / pb  is the quadratic sum of the
    # relative errors on psb and pb. Note that special care is required
    # in case of p=0.
    # relative error on p_sb:
    r_sb = delta_psb/psb if psb > 0 else 0
    r_b = delta_pb/pb if pb > 0 else 0
    delta_cls = cls * sqrt(r_sb**2 + r_b**2)
    return cls, delta_cls
    
def scan_s0_clsvalue_toys(b0, delta_b, nobs, s0min, s0max, nscan=10, ntoys = 1000):
    delta_s = (s0max - s0min) / nscan
    dclsvals, clsvals, svals = [], [], []
    for i in range(nscan): # i=0...nscan-1
        s0 = s0min + i*delta_s
        cls, delta_cls = get_cls_toys(s0, b0, delta_b, nobs, ntoys)
        print "CLs for s0=%.2f: %.3g+-%.3g" % (s0, cls, delta_cls)
        svals.append(s0)
        clsvals.append(cls)
        dclsvals.append(delta_cls)
    plot_xye(svals, clsvals, dclsvals, 'cls-vs-s-%.2f-%s-toys.pdf' % (delta_b, nobs))

    
def ex2d_i():
    #scan_s0_clsvalue_toys(5.2, 0.0, 6, 3., 10.)
    # one can also scan more precisely at the point where one suspects the limit to be:
    scan_s0_clsvalue_toys(5.2, 0.0, 6, 7., 7.4, 10, 10000)
    # "answer": seems to be consistent with result from 2.c. that the limit is at 7.2
    

#ex2d_i()
    
def ex2d_ii():
    # by iteratively adapting, one can find a range in which the limit is:
    #scan_s0_clsvalue_toys(5.2, 2.6, 6, 3., 10., 10, 10000)
    scan_s0_clsvalue_toys(5.2, 2.6, 6, 8., 9.0, 10, 10000)
    #
    # answer: the 95% C.L. limit is at about s=8.4, although even with ntoy=10000, the uncertainty
    # is large.
    
#ex2d_ii()




# 2.e.
# for i., commenting out the code in ex2.py leads to a limit of about 1.1--1.2.

# ii:
def ex2e_ii():
    model = build_shape_model(b_rate_unc = 0.1, b_shape_unc = True)
    limit, limit_error = get_observed_cls_limit(model)
    print "limit for mu with b rate/shape uncertainty: %.3f +- %.3f" % (limit, limit_error)
    
#ex2e_ii()


