'''
for general use for running an MCMC. to use:
from mcmc_class import MCMC

class MyProblem(MCMC)

 def __init__(self):
   super(MCMC, self).__init__()

that line will provide the following defaults:
   
self.nwalkers = 20
        self.nsteps = 2000
        self.rejection_fraction = 0.2
        self.burninsteps = 1000
        self.debug = debug
        self.threads = 1

must define guess, setup_data, and OTHER THINGS !!! DEBUG JCZ 250418

'''

from __future__ import print_function
import pandas as pd
from pudb import set_trace as pause
import emcee
import numpy as np
import corner
class MCMC(object):
    
    def __init__(self, *args, **kwargs):
        debug=kwargs.pop('debug', False)
        self.plot_dir = kwargs.pop('plot_dir', 'plots/')
        if debug:
            # JCZ 060818
            # was 10
            self.nwalkers = 20
            # JCZ 030518
            # !!! should be 100
            # !!! for purposes of DR2, should be 1000
            self.nsteps = 100
            self.rejection_fraction = 0.2
            self.burninsteps = 10
            self.threads = 1
        else:
            # JCZ 061220
            # should be 20 instead of 40
            self.nwalkers = 30
            # JCZ 310718
            # made to be 4000 instead of 100000
            # ran the table JCZ 301220 with 5000 and 30 walkers
            self.nsteps = 5000
            self.rejection_fraction = 0.2
            self.burninsteps = 1000
            self.threads = 1

    @property
    def names(self, *args, **kwargs):
        # JCZ 201117 
        if self._name_set: 
            return self._names 
        else:
            return ['']

    @names.setter
    def names(self, value):
        self._name_set = True
        self._names = value

    def setup_run(self, **kwargs):
        '''
        prepares the run. this is called within run().
        '''
        self.setup_data()
        # JCZ 200819
        # took out an a=20.0 term.
        self.mcmc = emcee.EnsembleSampler(self.nwalkers, self.nvariables, self.logl, threads=self.threads, kwargs=kwargs)
        
        self.start = ([tuple((self.guess)) + 1e-4*np.random.randn(self.nvariables) for i in range(self.nwalkers)])

    def run(self, **kwargs):
        '''
        kwargs are keyword arguments for self.logl
        '''
        self.setup_run(**kwargs)
        pos, prob, start = self.mcmc.run_mcmc(self.start, self.burninsteps)
        self.mcmc.reset()
        pos, prob, state = self.mcmc.run_mcmc(pos, self.nsteps)

        # defines the collapsed array of the parameter chain
        # self.samples
        # and their probabilities
        # self.probs
        self.wrapup_run(**kwargs)
    def save(self):
        # write the result to a file, based on the problem name.
        df_result = pd.DataFrame({})
        for val, err, name in zip(self.best_fit, self.error, self.names):
            df_result[name.replace('$', '')] = [val]
            df_result[name.replace('$', '')+'_err'] = [err]
        print(df_result)
        df_result.to_csv(self.problem + '_mcmc_best_fit.dat', index=False, sep=',')
        print( 'saved results to file {}'.format(self.problem + '_mcmc_best_fit.dat'))

    def wrapup_run(self, **kwargs):
        '''
        Inputs
        use_mean : bool
         if True, compute 
        '''
        use_mean = kwargs.get('use_mean', False)
        st = int(self.rejection_fraction*self.nsteps)
        end = -1

        self.samples = (self.mcmc.chain[:,st:end, :].reshape((-1, self.nvariables)))
        self.probs = (self.mcmc.lnprobability[:,st:end].reshape((-1)))
        # only calculate final values from those chains that don't have disallowed values
        good_ind = np.where(~np.isinf(np.abs(self.probs)))[0]
        # JCZ 170818
        # i think i need to also account for times when logl is nan... i've never considered that possibility before...
        good_ind = np.where(np.logical_and(~np.isinf(np.abs(self.probs)), ~np.isnan(np.abs(self.probs))))[0]
        self.samples = self.samples[good_ind,:]
        print( good_ind)
        print (self.samples)
        print( np.isnan(self.samples).any())
        print( np.isinf(self.samples).any())
        print( np.min(self.samples))
        print (np.min(self.samples[:,-1]))
        print (self.probs[np.argmin(self.samples[:,-1])])
        print( '^^^^')
        if len(self.samples) > 1:


            avg = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(self.samples, [16, 50, 84], axis=0)))

            # JCZ 010618
            # the two v definitions in the casesbelow used to be swtiched... that is fixed now.
            # also now both sets of definitions are saved.
            self.mean = np.mean(self.samples, axis=0)
            # with sqrt(N-1) in den.
            self.err_std = np.std(self.samples, axis=0 ,ddof=1.0)

            
            self.err_quant = np.array([a[1]/2. + a[2]/2. for a in avg])
            self.median = np.array([a[0] for a in avg])
            # JCZ 010121
            # trying out this definition because i don't think median is being defined
            self.median = np.median(self.samples, axis=0)
        else:
            raise Exception('there are 1 or fewer chains in the MCMC chain. retry.')









        if use_mean:
            self.best_fit = np.array(self.mean)
        else:
            self.best_fit = np.array(self.median)
        print ('set best_fit to :')
        print( self.best_fit)
        try:
            if use_mean:
                self.error = self.err_std.copy()
            else:
                self.error = self.err_quant.copy()
            print ('set error to:')
            print (self.error)
        except:
            print( 'cant set self.error --- maybe because it is a property of the parent class already.')

        self.acc_frac = np.mean(self.mcmc.acceptance_fraction)
        self.logl_max = np.max(self.mcmc.lnprobability)
        # self.max_fit = max_chain(self.mcmc.lnprobability, self.mcmc.chain, np.array(range(self.nvariables)))

        # print 'plotting the marginal distributions'
        # p_marg = Plot(figsize=(50,14),fname=outfile + '_dnu_marg', fdirect=MainParams.fdirect, mail=True, save=True)
        # p_marg.fig.clf()

        # for i in xrange(p_mcmc.nvariables):
        #     try:
        #         ax = p_marg.fig.add_subplot(p_mcmc.nvariables//3 + (p_mcmc.nvariables % 3 > 0), 3, i+1)
        #         ax.hist(samples[:,i], bins=int(np.sqrt(p_mcmc.nsteps - st)), normed=True)
        #         ax.vlines(output['best'][i] - output['error'][i], ax.get_ylim()[0], ax.get_ylim()[1], color='white', linestyle='dashed')
        #         ax.vlines(output['best'][i] + output['error'][i], ax.get_ylim()[0], ax.get_ylim()[1], color='white', linestyle='dashed')
        #         ax.set_ylabel(output['names'][i])
        #         ax.tick_params(labelsize=6)
        #         ax.tick_params(labelsize=6)
        #     except:
        #         print 'couldnt plot {}th parameter'.format(i)
        # p_marg.fig.tight_layout()
        # p_marg.wrapup()


        try:
            self.autocorr = np.array([a for a in self.mcmc.get_autocorr_time()])
        except:
            print ('could not create an autocorrelation time array because the chain was not long enough. setting autocorrelations for all variables to zero.')

            
        self.done = True
        
    def plot(self, show=False, save=True, ext='pdf'):
        '''
        plot a corner plot of the mcmc run.
        '''
        if not self.done:
            raise Exception('must run mcmc before trying to plot. do cluster.run().')
        
        print( np.min(self.samples[:,-1])        )

        # fig = corner.corner(self.samples, truths=self.best_fit,labels=self.names, label_kwargs={'pad':2.3, 'fontsize':35, 'labelsize':35}, plot_contours=True, use_math_text=True)
        # JCZ 280722
        # with new matplotlib have to change this i think:
        fig = corner.corner(self.samples, truths=self.best_fit,labels=self.names, label_kwargs={'fontsize':35}, plot_contours=True, use_math_text=True)
        ax = fig.gca()


        if save:
            

            
            
            fig.set_size_inches(19, 19)

            fig.subplots_adjust(bottom=0.2)
            fig.subplots_adjust(left=0.2)

            fig.savefig(self.plot_dir + self.problem+'_mcmc.'+ext, format=ext)
        if show:
            plt.show(); plt.ion()
