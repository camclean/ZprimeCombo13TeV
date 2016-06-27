# this example shows how to plot distributions for the fitted signal strength, its
# uncertainty and the "pull"

model = test_model.simple_counting(s = 200.0, b = 10000.0)
result = mle(model, input = 'toys:1', n= 100)

bs = []
delta_bs = []
pulls = []

for b, db in result['s']['beta_signal']:
    bs.append(b)
    delta_bs.append(db)
    pulls.append((1 - b)/db)
    
pdbs = plotdata()
pdbs.histogram(bs, 0.0, 2.0, 30, include_uoflow = True)
plot(pdbs, 'bs', 'ntoys', 'beta_signal.pdf')

pdd = plotdata()
pdd.histogram(delta_bs, 0.0, 1.0, 30, include_uoflow = True)
plot(pdd, 'dbs', 'ntoys', 'delta_beta_signal.pdf')

pdp = plotdata()
pdp.histogram(pulls, -3.0, 3.0, 100, include_uoflow = True)
plot(pdp, 'pull', 'ntoys', 'pull.pdf')
    

# to write the data to a file, use e.g.:
pdp.write_txt('pull.txt')

# to write it to a root file:
write_histograms_to_rootfile({'pull': pdp.histo(), 'bs': pdbs.histo(), 'delta_bs': pdd.histo()}, 'mle_histos.root')
    

