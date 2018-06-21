from montepython.likelihood_class import Likelihood
import os
import numpy as np
from scipy import interpolate
from lmfit import Minimizer, Parameters, report_fit
from scipy.linalg import block_diag
import pprint, pickle
import matplotlib.pyplot as plt
import time

#Lyman alpha likelihood by Maria Archidiacono and Riccardo Murgia

class Lya(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(data, {'output': 'mPk'})
        self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': 1.5*self.khmax}) 
        
        self.kh_0 = np.linspace(self.khmin, self.khmax, num=self.k_size)

        self.DATASET = "mike-hires" 
        self.zp = 4.5

        lcdm_points = 33    #number of grid points for the lcdm case (i.e. alpha=0, regardless of beta and gamma values)
        params_numbers = 3  #number of non-astro params (i.e. alpha,beta and gamma)

        file_path = os.path.join(self.data_directory, self.grid_file)
        if os.path.exists(file_path):
           with open(file_path, 'r') as grid_file:
                line = grid_file.readline()
                while line.find('#') != -1:
                    line = grid_file.readline()
                while (line.find('\n') != -1 and len(line) == 1):
                    line = grid_file.readline()
                for index in xrange(self.grid_size):
                    self.alphas[index] = float(line.split()[0])
                    self.betas[index] = float(line.split()[1])
                    self.gammas[index] = float(line.split()[2])
                    line = grid_file.readline()
        else
           raise io_mp.ConfigurationError('Error: grid file is missing')
           exit()

        X_real = np.zeros((len(alphas), params_numbers)) #real params

        for k in range(len(alphas)):  #real params
           X_real[k][0] = k12(alphas[k], betas[k], gammas[k])    #k_1/2
           X_real[k][1] = betas[k]
           X_real[k][2] = gammas[k]

        #for the normalization (see Alex's notes) #???
        self.a_min = min(X_real[:,0])
        self.b_min = min(X_real[:,1])
        self.g_min = min(X_real[:,2])
        self.a_max = max(X_real[:,0])
        self.b_max = max(X_real[:,1])
        self.g_max = max(X_real[:,2])

        #redshift independent parameters - params order: z_reio, sigma_8, n_eff, f_UV
        zind_param_size = [3, 5, 5, 3] #how many values I have for each param
        zind_param_min = np.array([7., 0.5, -2.6, 0.])
        zind_param_max = np.array([15., 1.5, -2.0, 1.])
        zind_param_ref = np.array([9., 0.829, -2.3074, 0.])
        zreio_range = zind_param_max[0]-zind_param_min[0]
        neff_range = zind_param_max[2]-zind_param_min[2]

        # redshift dependent parameters - params order: params order: mean_f , t0, slope
        zdep_params_size = [9, 3, 3] #how many values I have for each param
        zdep_params_refpos = [4, 1, 2] #where to store the P_F(ref)# DATA 

        #MEAN FLUXES values###
        flux_ref_old = (np.array([0.669181, 0.617042, 0.564612, 0.512514, 0.461362, 0.411733, 0.364155, 0.253828, 0.146033, 0.0712724]))
        flux_min_meanf = (np.array([0.401509, 0.370225, 0.338767, 0.307509, 0.276817, 0.24704, 0.218493, 0.152297, 0.0876197, 0.0427634]))
        flux_max_meanf = (np.array([0.936854, 0.863859, 0.790456, 0.71752, 0.645907, 0.576426, 0.509816, 0.355359, 0.204446, 0.0997813]))

        # DATA 
        # FIRST (NOT USED) DATASET (19 wavenumbers) ***XQ-100***
        zeta_range_XQ = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2]  #list of redshifts corresponding to the 19 wavenumbers (k) (see ./SPECTRA/README)
        k_XQ = [0.003,0.006,0.009,0.012,0.015,0.018,0.021,0.024,0.027,0.03,0.033,0.036,0.039,0.042,0.045,0.048,0.051,0.054,0.057]

        # SECOND DATASET (7 wavenumbers) ***HIRES/MIKE***
        zeta_range_hm = [4.2, 4.6, 5.0, 5.4]  #list of redshifts corresponding to the 7 wavenumbers (k) (see ./SPECTRA/README)
        k_hm = [0.00501187,0.00794328,0.0125893,0.0199526,0.0316228,0.0501187,0.0794328] #in s/km

        self.zeta_full_lenght = (len(zeta_range_XQ) + len(zeta_range_hm))
        self.kappa_full_lenght = (len(k_XQ) + len(k_hm))
        redshift = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.2, 4.6, 5.0, 5.4] #which snapshots (first 7 for first dataset, last 4 for second one)

        #T0 AND SLOPE VALUES
        t0_ref_old = np.array([11251.5, 11293.6, 11229.0, 10944.6, 10421.8, 9934.49, 9227.31, 8270.68, 7890.68, 7959.4])
        slope_ref_old = np.array([1.53919, 1.52894, 1.51756, 1.50382, 1.48922, 1.47706, 1.46909, 1.48025, 1.50814, 1.52578])

        t0_values_old = np.zeros(( 10, zdep_params_size[1] ))
        t0_values_old[:,0] = np.array([7522.4, 7512.0, 7428.1, 7193.32, 6815.25, 6480.96, 6029.94, 5501.17, 5343.59, 5423.34])
        t0_values_old[:,1] = t0_ref_old[:]
        t0_values_old[:,2] = np.array([14990.1, 15089.6, 15063.4, 14759.3, 14136.3, 13526.2, 12581.2, 11164.9, 10479.4, 10462.6])

        slope_values_old = np.zeros(( 10, zdep_params_size[2] ))
        slope_values_old[:,0] = np.array([0.996715, 0.979594, 0.960804, 0.938975, 0.915208, 0.89345, 0.877893, 0.8884, 0.937664, 0.970259])
        slope_values_old[:,1] = [1.32706, 1.31447, 1.30014, 1.28335, 1.26545, 1.24965, 1.2392, 1.25092, 1.28657, 1.30854]
        slope_values_old[:,2] = slope_ref_old[:]

        t0_min = t0_values_old[:,0]*0.1; t0_max = t0_values_old[:,2]*1.4
        slope_min = slope_values_old[:,0]*0.8; slope_max = slope_values_old[:,2]*1.15

        #IMPORTING THE TWO GRIDS FOR KRIGING
        # Here I import the gridS that I pre-computed through the file "setting_Kriging_grid_2R.py"  
        file_path = os.path.join(self.data_directory, self.astro_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           full_matrix_interpolated_ASTRO = pickle.load(pkl)
           print full_matrix_interpolated_ASTRO.shape
        else
           raise io_mp.ConfigurationError('Error: astro spectra file is missing')
           exit()

        file_path = os.path.join(self.data_directory, self.abg_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           self.full_matrix_interpolated_ABG = pickle.load(pkl)
           print self.full_matrix_interpolated_ABG.shape
        else
           raise io_mp.ConfigurationError('Error: abg spectra file is missing')
           exit()

        ALL_zdep_params = len(flux_ref_old) + len(t0_ref_old) + len(slope_ref_old)
        grid_lenght_ABG = len(self.full_matrix_interpolated_ABG[0,0,:])
        grid_lenght_ASTRO = len(full_matrix_interpolated_ASTRO[0,0,:])
        astroparams_number_KRIG = len(zind_param_size) + ALL_zdep_params

        #### --- ABG GRID --- {alpha, beta, gamma} GRID
        file_path = os.path.join(self.data_directory, self.abg_grid_file)
        if os.path.exists(file_path):
           self.X_ABG = np.zeros((grid_lenght_ABG, params_numbers))
           for param_index in range(params_numbers):
               self.X_ABG[:,param_index] = np.genfromtxt(file_path, usecols=[param_index], skip_header=1)
        else
           raise io_mp.ConfigurationError('Error: abg grid file is missing')

        #### --- ASTRO GRID --- ORDER OF THE COMPLETE LIST OF ASTRO PARAMS: z_reio, sigma_8, n_eff, f_UV, mean_f(z), t0(z), slope(z)
        #### HIGH REDSHIFT
        file_path = os.path.join(self.data_directory, self.abg+astro_grid_file)
        if os.path.exists(file_path):
           X = np.zeros((grid_lenght_ASTRO,astroparams_number_KRIG))
           for param_index in range(astroparams_number_KRIG):
               X[:,param_index] = np.genfromtxt(file_path, usecols=[param_index], skip_header=1)
        else
           raise io_mp.ConfigurationError('Error: abg+astro grid file is missing')

        # FUNCTIONS FOR THE ORDINARY KRIGING INTERPOLATION ON ALL THE (k,z) COMBINATIONS #################################
        self.epsilon = 1e-8
        self.exponent = 6.
        #  STUFF FOR INTERPOLATING IN THE ASTROPARAMS SPACE ####################################
        redshift_list = np.array([3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.4]) #combined dataset (MIKE/HIRES + XQ-100)
        astrokrig_result = np.zeros((self.zeta_full_lenght, self.kappa_full_lenght ))
        #minimum and maximum values for the kriging normalisation###
        self.F_prior_min = np.array([0.535345,0.493634,0.44921,0.392273,0.338578,0.28871,0.218493,0.146675,0.0676442,0.0247793])
        self.F_prior_max = np.array([0.803017,0.748495,0.709659,0.669613,0.628673,0.587177,0.545471,0.439262,0.315261,0.204999])

        #DATA
        model_H = np.zeros (( len(zeta_range_mh), len(k_mh) ))
        y_H = np.zeros (( len(zeta_range_mh), len(k_mh) ))
        model_M = np.zeros (( len(zeta_range_mh)-1, len(k_mh) ))
        y_M = np.zeros (( len(zeta_range_mh)-1, len(k_mh) ))

        file_path = os.path.join(self.data_directory, self.MIKE_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           y_M_reshaped = pickle.load(pkl)
        else
           raise io_mp.ConfigurationError('Error: MIKE spectra file is missing')
        
        file_path = os.path.join(self.data_directory, self.HIRES_spectra_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           y_H_reshaped = pickle.load(pkl)
        else
           raise io_mp.ConfigurationError('Error: HIRES spectra file is missing')

        file_path = os.path.join(self.data_directory, self.MIKE_cov_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           cov_M_inverted = pickle.load(pkl)
        else
           raise io_mp.ConfigurationError('Error: MIKE covariance matrix file is missing')

        file_path = os.path.join(self.data_directory, self.HIRES_cov_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           cov_H_inverted = pickle.load(pkl)
        else
           raise io_mp.ConfigurationError('Error: HIRES covariance matrix file is missing')

        file_path = os.path.join(self.data_directory, self.PF_noPRACE_file)
        if os.path.exists(file_path):
           pkl = open(file_path, 'r')
           self.PF_noPRACE = pickle.load(pkl)
        else
           raise io_mp.ConfigurationError('Error: PF_noPRACE file is missing')

        cov_MH_inverted = block_diag(cov_H_inverted,cov_M_inverted)
        y_MH_reshaped = np.concatenate((y_H_reshaped, y_M_reshaped))

        return



    #from alpha to 1./k_{1/2} in order to interpolate in a less sparse grid
    def k12(alpha,beta,gamma):
        return ((((0.5)**(1/(2*gamma)) - 1)**(1/beta))/alpha)**(-1)       #1./k1half



    def loglkl(self, cosmo, data):

        #deal with the astro nuisance parameters
        if 'T0a' in self.use_nuisance:
            T0a=data.mcmc_parameters['T0a']['current']*data.mcmc_parameters['T0a']['scale']
        if 'T0s' in self.use_nuisance:
            T0s=data.mcmc_parameters['T0s']['current']*data.mcmc_parameters['T0s']['scale']
        if 'gamma_a' in self.use_nuisance:
            gamma_a=data.mcmc_parameters['gamma_a']['current']*data.mcmc_parameters['gamma_a']['scale']
        if 'gamma_s' in self.use_nuisance:
            T0s=data.mcmc_parameters['gamma_s']['current']*data.mcmc_parameters['gamma_s']['scale']
        if 'Fz1' in self.use_nuisance:
            Fz_names=[]
            for index_z in xrange(4):
                Fz_names.append('Fz'+str(index_z+1))
                Fz[index_z+1] = data.mcmc_parameters[Fz_names[index_z]]['current']*data.mcmc_parameters[Fz_names[index_z]]['scale']
        if 'F_UV' in self.use_nuisance:
            F_UV=data.mcmc_parameters['F_UV']['current']*data.mcmc_parameters['F_UV']['scale']

        h=data.cosmo_arguments['H0']/100.
        Plin_0 = np.zeros(len(self.kh_0), 'float64')
        for index_k in range(len(self.kh_0)):
            Plin_0[index_k] = cosmo.pk(self.kh_0[index_k]*h, 0.0) #use pk_lin with the new class version 
        Plin_0 *= h**3


        #see likelihood_class get_flat_fid
        param_backup = data.cosmo_arguments
        print '\n'
        print 'param_backup'
        print param_backup
        #pba->Omega0_idr = pba->stat_f_idr*pow(pba->xi_idr,4.)*pba->Omega0_g;
        #N_dark = pba->Omega0_idr/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;
        DeltaNeff=data.cosmo_arguments['stat_f_idr']*(data.cosmo_arguments['xi_idr']**4)/7.*8./((4./11.)**(4./3.))
        eta2=(1.+0.2271*(data.cosmo_arguments['N_ur']+DeltaNeff))/(1.+0.2271*data.cosmo_arguments['N_ur'])
        eta=np.sqrt(eta2)
        ob_backup = data.cosmo_arguments['omega_b']
        data.cosmo_arguments['omega_b'] *= 1./eta2
        oc_backup = data.cosmo_arguments['omega_cdm']
        data.cosmo_arguments['omega_cdm'] *= 1./eta2
        H0_backup = data.cosmo_arguments['H0']
        data.cosmo_arguments['H0'] *= 1./eta
        h=data.cosmo_arguments['H0']/100.
        xi_backup = data.cosmo_arguments['xi_idr']
        data.cosmo_arguments['xi_idr'] = 0.
        f_backup = data.cosmo_arguments['f_idm_dr']
        data.cosmo_arguments['f_idm_dr'] = 0.
        a_backup = data.cosmo_arguments['a_dark']
        data.cosmo_arguments['a_dark'] = 0.
        print 'DeltaNeff = ',DeltaNeff,' eta^2 = ',eta2
        print '\n'
        #convert Neff according to arXiv:
        #data.cosmo_arguments = {'output': ' mPk ','P_k_max_h/Mpc': 1.5*self.khmax,
        #                        'omega_b': ob,'omega_cdm': oc,'H0':H0,
        #                        'A_s': data.cosmo_arguments['A_s'],'n_s': data.cosmo_arguments['n_s'], 'tau_reio': data.cosmo_arguments['tau_reio']        #                        ,'N_ur':'3.046'}

        cosmo.empty()
        print 'equivalent lcdm'
        print data.cosmo_arguments
        cosmo.set(data.cosmo_arguments)
        cosmo.compute(['lensing'])

        Plin_equiv_0 = np.zeros(len(self.kh_0), 'float64')
        for index_k in range(len(self.kh_0)):
            Plin_equiv_0[index_k] = cosmo.pk(self.kh_0[index_k]*h, 0.0) #use pk_lin with the new class version
        Plin_equiv_0 *= h**3
        z_reio=cosmo.z_reio()
        neff=cosmo.neff()
        sigma8=cosmo.sigma8()
        print 'z_reio = ',z_reio,'sigma8 = ',sigma8,' neff = ',neff
        print '\n'

        cosmo.empty()
        data.cosmo_arguments = param_backup #this does not work
        data.cosmo_arguments['omega_b'] = ob_backup 
        data.cosmo_arguments['omega_cdm'] = oc_backup
        data.cosmo_arguments['H0'] = H0_backup
        h=data.cosmo_arguments['H0']/100.
        data.cosmo_arguments['xi_idr'] = xi_backup
        data.cosmo_arguments['f_idm_dr'] = f_backup
        data.cosmo_arguments['a_dark'] = a_backup
        print 'back to the original model'
        print data.cosmo_arguments
        print '\n'
        cosmo.set(data.cosmo_arguments)
        cosmo.compute(['lensing'])

        Tk_0 = np.zeros(len(self.kh_0), 'float64')
        Tk_0 = np.sqrt(Plin_0/Plin_equiv_0) 
        spline=interpolate.splrep(self.kh_0,Tk_0)
        der = interpolate.splev(self.kh_0, spline, der=1)

        #Now merge with Riccardo's interpolation code
        #setting k_max (i.e. cutting oscillations from the fitted region)
        for index_k in range(len(self.kh_0)):
            index_khmax = -1
            if Tk_0[index_k]<0.2 and der[index_k]>0.:
               index_khmax = index_k
               print index_khmax
               break

        kh = self.kh_0[:index_khmax]
        Plin_equiv = Plin_equiv_0[:index_khmax]
        Plin = Plin_0[:index_khmax]

        Tk = np.sqrt(Plin/Plin_equiv)

        # fitting the given linear P(k) with the {alpha,beta,gamma}-formula 

        #model function #T^2=P_model/P_ref
        def T(kh,alpha,beta,gamma):
            return (1. + (alpha*kh)**(beta))**(gamma)
        #define objective function: returns the array to be minimized
        def fcn2min(params, kh, Tk):
            alpha = params['alpha']
            beta = params['beta']
            gamma = params['gamma']
            model = T(kh,alpha,beta,gamma) #(1. + (alpha*kappa_interp)**(beta))**(gamma)
            return (model - Tk)      #standard residuals

        # create a set of Parameters
        params = Parameters()
        params.add('alpha', value=0.001, min = 0., max = 0.3)
        params.add('beta', value=2.24, min = 0.5, max = 10.)
        params.add('gamma', value=-4.46, min=-10., max=-0.1)

        # do fit, default is with least squares method
        t0_fit = time.clock()

        minner = Minimizer(fcn2min, params, fcn_args=(kh, Tk))
        result = minner.minimize(method = 'leastsq')
        best_alpha = result.params['alpha'].value
        best_beta  = result.params['beta'].value
        best_gamma = result.params['gamma'].value

        t1_fit = time.clock()

        # write error report
        report_fit(result)

        plt.xlabel('k [h/Mpc]')
        plt.ylabel('$P_{nCDM}/P_{CDM}$')

        plt.ylim(0.,1.1)
        plt.xlim(self.khmin,self.khmax)
        plt.xscale('log')
        #plt.yscale('log')
        plt.grid(True)

        plt.plot(self.kh_0, Tk_0**2, 'r')
        plt.plot(self.kh_0, (T(self.kh_0, best_alpha, best_beta, best_gamma))**2, 'b--')
        #plt.plot(self.kh_0, der, 'k')
        #plt.show()
        plt.savefig('grid_fit_plot.pdf')

        ABG_matrix_new = np.zeros((self.zeta_full_lenght, self.kappa_full_lenght, len(self.alphas)+8 )) #??? why + 8? what is 8? define a parameter/quantity related to this 8

        NEW_ABG_matrix = np.zeros(( len(self.alphas)+8, self.zeta_full_lenght, self.kappa_full_lenght))
        for i in range(zeta_full_lenght):
            for j in range(kappa_full_lenght):
                NEW_ABG_matrix[:,i,j] = self.full_matrix_interpolated_ABG[i,j,:]

        def ordkrig_estimator_3D(p21, z):
            ABG_matrix_new = NEW_ABG_matrix + ordkrig_estimator(p21,z) - 1.
            ABG_matrix_new = np.clip(ABG_matrix_new, 0. , None)
            for i in range(self.zeta_full_lenght):
                for j in range(self.kappa_full_lenght):
                    self.full_matrix_interpolated_ABG[i,j,:] = ABG_matrix_new[:,i,j]
        return np.sum(np.multiply(ordkrig_lambda_3D((k12(p21[0],p21[1],p21[2]))/(self.a_max-self.a_min), p21[1]/(self.b_max-self.b_min), p21[2]/(self.g_max-self.g_min), self.X_ABG[:,0], self.X_ABG[:,1], self.X_ABG[:,2]), self.full_matrix_interpolated_ABG[:,:,:]),axis=2)

        def ordkrig_estimator(p21, z):
            pa10 = (z_dep_func(p21[11], p21[12], z[:])*1e4)/(t0_max[:]-t0_min[:])
            pb10 = z_dep_func(p21[13], p21[14], z[:])/(slope_max[:]-slope_min[:])
            p37 = np.concatenate((p21[3:11], pa10[6:], pb10[6:]))
            for index in range(7,len(redshift)): #??? 7 what is 7?
                astrokrig_result[index,:] = np.sum(np.multiply(ordkrig_lambda(p37[0]/zreio_range, p37[1], p37[2]/neff_range, p37[3], p37[4+index-7]/(F_prior_max[index-1]-F_prior_min[index-1]), p37[8+index-7], p37[12+index-7], X[:,0], X[:,1], X[:,2], X[:,3], X[:,4+index-1], X[:,14+index-1], X[:,24+index-1]), full_matrix_interpolated_ASTRO[index,:,:]),axis=1)
        return astrokrig_result

        def z_dep_func(parA, parS, z):  #analytical function for the redshift dependence of t0 and slope
        return parA*(( (1.+z)/(1.+self.zp) )**parS)


        model = self.PF_noPRACE*ordkrig_estimator_3D(theta, z)
        upper_block = np.vsplit(model, [7,11])[0]
        lower_block = np.vsplit(model, [7,11])[1]
        if self.DATASET == "mike-hires":
           model_H[:,:] = lower_block[:,19:]; model_H_reshaped = np.reshape(model_H, -1, order='C')
           model_M[:,:] = lower_block[:3,19:]; model_M_reshaped = np.reshape(model_M, -1, order='C')
           model_MH_reshaped = np.concatenate((model_H_reshaped,model_M_reshaped))
           chi2 = np.dot((y_MH_reshaped - model_MH_reshaped),np.dot(cov_MH_inverted,(y_MH_reshaped - model_MH_reshaped)))
        else
           raise io_mp.LikelihoodError('Error: for the time being, only the mike - hires dataset is available')
           exit()

        print 'chi^2 = ',chi2
        return -chi2/2.
