model = test_model.simple_counting(s = 1.0, n_obs = 10.0, b = 5.0, b_uncertainty = 3.0)
model_summary(model)
result = mle(model, input = 'data', n = 1)
print "result: ", result
report.write_html('htmlout')
