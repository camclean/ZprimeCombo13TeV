#!/usr/bin/env python

channels = ['muo', 'ele']
signals  = ['zpn', 'zpw', 'zph', 'kkg']

sig_mass = {
  'zpn': [500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000],
  'zpw': [500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000],
  'zph': [500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000],
  'kkg': [500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000],
}
####

onlyCR = 0
onlySR = 1

def obs_veto(h_):
    veto = False
#    if 'mtt' not in h_.split('__')[0]: veto = True

    if onlyCR and 'L1chi2lo'     in h_.split('__')[0]: veto = True
    if onlySR and 'L1chi2lo' not in h_.split('__')[0]: veto = True

#    if 't0b0'     in h_.split('__')[0]: veto = True
#    if 't0B1'     in h_.split('__')[0]: veto = True
#    if 'T1B0'     in h_.split('__')[0]: veto = True

#    if 'muo' not in h_.split('__')[0]: veto = True
#    if 'ele' not in h_.split('__')[0]: veto = True

    return veto

def analysis_def(key):
    if   key == 'SIG'    : return [[s+'*' for s in signals], SIG]
    elif key == 'SIG_muo': return [[s+'*' for s in signals], SIG_muo]
    elif key == 'SIG_ele': return [[s+'*' for s in signals], SIG_ele]

    elif key == 'zpn'    : return ['zpn*', ZPN]
    elif key == 'zpn_muo': return ['zpn*', ZPN_muo]
    elif key == 'zpn_ele': return ['zpn*', ZPN_ele]

    elif key == 'zpw'    : return ['zpw*', ZPW]
    elif key == 'zpw_muo': return ['zpw*', ZPW_muo]
    elif key == 'zpw_ele': return ['zpw*', ZPW_ele]

    elif key == 'zph'    : return ['zph*', ZPH]
    elif key == 'zph_muo': return ['zph*', ZPH_muo]
    elif key == 'zph_ele': return ['zph*', ZPH_ele]

    elif key == 'kkg'    : return ['kkg*', KKG]
    elif key == 'kkg_muo': return ['kkg*', KKG_muo]
    elif key == 'kkg_ele': return ['kkg*', KKG_ele]

    else:
        print '\n @@@ FATAL -- undefined model key: '+key
        raise SystemExit

# models
def SIG    (h): return (True                        and (not obs_veto(h)))
def SIG_muo(h): return (True             and muo(h) and (not obs_veto(h)))
def SIG_ele(h): return (True             and ele(h) and (not obs_veto(h)))

def ZPN    (h): return (signal('zpn', h)            and (not obs_veto(h)))
def ZPN_muo(h): return (signal('zpn', h) and muo(h) and (not obs_veto(h)))
def ZPN_ele(h): return (signal('zpn', h) and ele(h) and (not obs_veto(h)))

def ZPW    (h): return (signal('zpw', h)            and (not obs_veto(h)))
def ZPW_muo(h): return (signal('zpw', h) and muo(h) and (not obs_veto(h)))
def ZPW_ele(h): return (signal('zpw', h) and ele(h) and (not obs_veto(h)))

def ZPH    (h): return (signal('zph', h)            and (not obs_veto(h)))
def ZPH_muo(h): return (signal('zph', h) and muo(h) and (not obs_veto(h)))
def ZPH_ele(h): return (signal('zph', h) and ele(h) and (not obs_veto(h)))

def KKG    (h): return (signal('kkg', h)            and (not obs_veto(h)))
def KKG_muo(h): return (signal('kkg', h) and muo(h) and (not obs_veto(h)))
def KKG_ele(h): return (signal('kkg', h) and ele(h) and (not obs_veto(h)))

def signal(sig, h_):
    if sig not in signals:
        print '\n @@@ FATAL -- undefined signal key: '+sig
        raise SystemExit

    prc = h_.split('__')[1]
    for s in filter(lambda a: a != sig, signals):
        if s in prc: return False

    if sig in prc:
        mass = prc.replace(sig, '')
        return float(mass) in sig_mass[sig]

    return True

# channels
def muo(h): return obs_prex(h, 'muo')
def ele(h): return obs_prex(h, 'ele')

def obs_prex(h_, obs_prex_):
    return h_.split('__')[0].startswith(obs_prex_)
####

def build_model__ttbar_ljets(files, hfilter, signal, mcstat):
    mod = build_model_from_rootfile(files, hfilter, include_mc_uncertainties = mcstat)
    mod.fill_histogram_zerobins()
    mod.set_signal_processes(signal)

    for p in mod.processes: mod.add_lognormal_uncertainty('lumi', math.log(1.027), p)

    mod.add_lognormal_uncertainty('xsec_ttbar', math.log(1.08), 'ttbar')
    mod.add_lognormal_uncertainty('xsec_wjets', math.log(1.06), 'wjetL')
    mod.add_lognormal_uncertainty('xsec_wjets', math.log(1.06), 'wjetH')
    mod.add_lognormal_uncertainty('xsec_sitop', math.log(1.20), 'sitop')
    mod.add_lognormal_uncertainty('xsec_zjets', math.log(1.20), 'zjets')
    mod.add_lognormal_uncertainty('xsec_dibos', math.log(1.20), 'dibos')

    return mod

def build_model(files_, key_, mcstat_ = True):

    sig, hfilter = analysis_def(key_)

    model = build_model__ttbar_ljets(files_, hfilter, sig, mcstat_)

    free_params = [
#      'xsec_ttbar',
#      'xsec_wjets',
#      'xsec_zjets',
      'ttagT',
    ]

    fixd_params = [
#      'muRF_ttbar',
#      'muRF_wjets',
    ]

    for p in model.distribution.get_parameters():
        if p in free_params: model.distribution.set_distribution_parameters(p, width = float('inf'))
        if p in fixd_params: model.distribution.set_distribution_parameters(p, width = float(.0001))

    return model

# theta__run add-on
ifiles = ["xttlj13__160218_T1_v01__l2p60__160222__mod__rebinned_e30_m50__rescaled.root"]
an_key = "zpn"
model = build_model(ifiles, an_key, True)

model_summary(model)

opts = Options()

RUN_THETA = True

if not RUN_THETA:
    exp_l = bayesian_quantiles(model, input='toys:0', n=1000, run_theta=False, hint_method='zero')
    obs_l = bayesian_quantiles(model, input='data'  , n=  10, run_theta=False, hint_method='zero')
    for spid in exp_l: exp_l[spid].get_configfile(opts)
    for spid in obs_l: obs_l[spid].get_configfile(opts)
else:
    res = bayesian_limits(model, what='all')
    exp, obs = res

    o_file = open('limits.txt', 'w')
    for i in range(len(exp.x)):
        o_file.write( '%.2f %.5f' % (exp.x[i], exp.y[i]))
        o_file.write(' %.5f %.5f' % (exp.bands[1][1][i], exp.bands[1][0][i]))
        o_file.write(' %.5f %.5f' % (exp.bands[0][1][i], exp.bands[0][0][i]))
        o_file.write(' %.5f'      % (obs.y[i] if obs else -1.))
        o_file.write('\n')
    o_file.close()
