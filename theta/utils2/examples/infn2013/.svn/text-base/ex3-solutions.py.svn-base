execfile("ex3.py")
    
def ex3a_i():
    svals, ps = counting_posterior(6, 5.2, 0.0, 10.0)
    plot_xy(svals, ps, 'counting_posterior.pdf')
    # answer: it is at 6-5.2 = 0.8 (can also be shown analytically for this example!)

    
#ex3a_i()


def ex3a_ii():
    svals, ps = counting_posterior_unc(6, 5.2, 2.6, 0.0, 10.0)
    plot_xy(svals, ps, 'counting_posterior_unc.pdf')
    # answer: the posterior is broader; the maximum posterior is at s=0 now, so it
    # has shifted.

#ex3a_ii()



# 3.b.

def ex3b():
    model = build_shape_model()
    mus, post = get_posterior(model)
    plot_xy(mus, post, 'shape_posterior0.pdf')
    model = build_shape_model(b_rate_unc = 0.1, b_shape_unc = True)
    mus, post = get_posterior(model)
    plot_xy(mus, post, 'shape_posterior_unc.pdf')
    # answer: maximum stays approximately the same, but width increases.

#ex3b()
    
    
   
# 3.c.    

def get95up(xpost, ypost):
    ytotal = sum(ypost)
    s = 0.0
    for i in range(len(xpost)):
        s += ypost[i]
        if s/ytotal > 0.95: return xpost[i]
    

def ex3c_i():
    x, y = counting_posterior(6, 5.2, 0.0, 20.0, 100)
    l = get95up(x, y)
    print "counting experiment b=5.2; nobs=6"
    print "Bayesian limit for s: %.3g" % l
    x, y = counting_posterior_unc(6, 5.2, 2.6, 0.0, 20.0, 100)
    l = get95up(x, y)
    print "counting experiment b=5.2 +- 2.6; nobs=6"
    print "Bayesian limit for s: %.3g" % l
    
    # result:
    # limits are 6.9 and 7.5, resp.
    
#ex3c_i()


def ex3c_ii():
    x, y = counting_posterior(0, 0, 0.0, 10.0, 1000)
    l = get95up(x, y)
    print "counting experiment b=0.0; nobs=0"
    print  "Bayesian limit for s: %.3g" % l
    plot_xy(x, y, 'post0_0.pdf')
    # answer: the limit is 3.0. It coincides with the frequentist limit.

    
#ex3c_ii()



# 3.d.

def ex3d():
    x, y = counting_posterior(6, 5.2, 0.0, 20.0, 1000)
    plot_xy(x, y, 'post6_52-flatprior.pdf')
    l = get95up(x, y)
    print  "Bayesian limit for flat prior: %.3g" % l
    apply_prior(x, y)
    l = get95up(x, y)
    plot_xy(x, y, 'post6_52-linearprior.pdf')
    print  "Bayesian limit for linear prior: %.3g" % l
    # answer: limit changes from 7.2 to 9.4
    
#ex3d()
    
