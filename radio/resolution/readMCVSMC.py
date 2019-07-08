import cPickle
import numpy
import matplotlib.pyplot as plt
import numpy as np


event='223619964'

methfile='meth223619964.dat'
f=open(methfile)
mcvsmc_info = cPickle.load(f) # contains tuple

reco_file=open('reco95102507.dat',"r")
reco_info=cPickle.load(reco_file)
reco_file.close
print reco_info.keys()

(nsimprot, nsimiron, ev, xreco, xreal, xbest, cchi2, rchi2, p_ratio, p_ratio0, p_ratio1, d_ratio, dratio, sim_tot_power, rbf, d150, nch, mr, sa, xoff, yoff, xcore, ycore, en, zen, az, noisepower, antratio, nsim, nsim_prot) = mcvsmc_info

selection_valid_reconstructions = (xreco > 1)
nfailed = np.sum(xreco < 1)

sel = selection_valid_reconstructions


xmax_reco = reco_info['xmaxreco']

#print d_ratio
#print '___________________'
#print p_ratio

#xmax_reco = reco_info['xmaxreco']
sigma_e = np.std(d_ratio[sel])
sigma_e_radio = np.std(np.sqrt(p_ratio[sel]))
sigma_logE_radio = np.std(0.5 * np.log10(p_ratio[sel]))


'''
fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)


bins=np.arange(0.2,1.8,0.05)

ax1.hist(p_ratio.flatten(),bins=bins,histtype='step',color='blue')
ax1.hist(p_ratio[sel].flatten(),bins=bins,histtype='step',color='green')

ax2.hist(d_ratio.flatten(),bins=bins,histtype='step',color='blue')
ax2.hist(d_ratio[sel].flatten(),bins=bins,histtype='step',color='green')


plt.show()
'''
'''
selection_valid_reconstructions = (xreco > 1)
nfailed = np.sum(xreco < 1)
sel = selection_valid_reconstructions

xmax_reco = reco_info['xmaxreco']
sigma_e = np.std(d_ratio[sel])
sigma_e_radio = np.std(np.sqrt(p_ratio[sel]))
sigma_logE_radio = np.std(0.5 * np.log10(p_ratio[sel]))

rms_err_x = np.sqrt( np.average( (xreal[sel] - xreco[sel]).ravel()**2 ))
'''

'''
p_ratio=info['p_ratio']
d_ratio=info['d_ratio']


d_ratios=info['d_ratios']
lofar_sigma=info['lofar_sigma']
pratio1=info['pratio1']
p_ratio0=info['p_ratio0']
print 'p_ratio: {0}'.format(p_ratio)

'''
'''
print 'p_ratio: {0}'.format(p_ratio)
print 'd_ratio: {0}'.format(d_ratio)

print 'd_ratios: {0}'.format(d_ratios)
print 'lofar_sigma: {0}'.format(lofar_sigma)
print 'pratio1: {0}'.format(pratio1)
print 'p_ratio0: {0}'.format(p_ratio0)
'''

