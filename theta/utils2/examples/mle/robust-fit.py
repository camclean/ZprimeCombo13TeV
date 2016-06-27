model = test_model.bernstein_model(10) # a model of 10 bernstein polynomials, difficult to minimize
mle(model, 'toys:0.0', 100)
mle(model, 'toys:0.0', 100, with_error = False)
options = Options()
options.set('minimizer', 'strategy', 'robust')
mle(model, 'toys:0.0', 100, with_error = False, options = options)

# note: in this particular case, the 'robust' strategy does not really help much, as with_error is enough

