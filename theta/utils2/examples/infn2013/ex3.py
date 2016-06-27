execfile("common.py")

# 3.a.

# returns x,y data (= s values, posterior density values)
# for the posterior density for s for a counting
# experiment with known background b and a flat prior for s for
# nobs observed events by scanning s from smin to smax, using nscan points.
#
# For the normalization, it is assumed that the posterior is negligible outside the
# range [smin, smax]
def counting_posterior(nobs, b, smin, smax, nscan=100):
    delta_s = (smax - smin) / nscan
    x, y = [], []
    for i in range(nscan):
        s = smin + delta_s * i
        x.append(s)
        # the posterior at this point is the poisson probability; note
        # that at this point, we don't care about normalization:
        y.append(poisson_p_eq(nobs, s + b))
    # now, scale the y values such that we have actually a density:
    factor = 1.0 / (sum(y) * delta_s)
    y = [factor * y0 for y0 in y]
    return x, y
        
#svals, ps = counting_posterior() # TODO: fill the blanks
#plot_xy(svals, ps, 'counting_posterior.pdf')



# return x,y data (= s vales, posterior density) for a model with
# an uncertainty on b, i.e. where b has a (truncated) normal prior
# with mean b0 and width delta_b.
# Note that the data for the posterior is NOT normalized.
def counting_posterior_unc(nobs, b0, delta_b, smin, smax, nscan = 100):
    # integrate over b in 100 steps from -5sigma to 5sigma:
    nscan_b = 100
    bmin, bmax = b0 - 5*delta_b, b0 + 5*delta_b
    # avoid 0 and negative values for b, as Poisson for negative means are not defined:
    if bmin <= 0: bmin = bmax / nscan_b # a small value on the scale we are interested in
    # the marginal posterior in s:
    ps = [0.0] * nscan
    # integration loop for b:
    for ib in range(nscan_b):
        b = bmin + ib * (bmax - bmin) / nscan_b
        # we use a normal prior for b. Note that we can neglect constant  factors
        # ("constant"=independent of b and s)
        b_prior = exp(-0.5 * (b - b0)**2 / delta_b**2)
        x, ps_b = counting_posterior(nobs, b, smin, smax, nscan)
        # calculate ps += b_prior * ps_b:
        for i in range(nscan):
            ps[i] += b_prior * ps_b[i]
    return x, ps


#svals, ps = counting_posterior_unc() # TODO: fill the blanks
#plot_xy(svals, ps, 'counting_posterior_unc.pdf')


# 3.b.

#model = build_shape_model()
#mus, post = get_posterior(model)
#plot_xy(mus, post, 'shape_posterior.pdf')



# 3.c.

# return the 95% C.L. upper limit, given the posterior as (x,y) values; x is a sorted list
# of equidistant values and y[i] contains the posterior density for this value of x[i].
def get95up(xpost, ypost):
    # The method should not assume that the y-values are normalized
    ytotal = sum(ypost)
    #for i in range(len(xpost)): TODO: implement!
        
    
    
# 3.d.

# updates the (x,y) posterior data by multiplying it with a prior p(x)=x
def apply_prior(x,y):
    for i in range(len(x)):
        y[i] *= x[i]




