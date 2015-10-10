import numpy as np
import h5py
import Analysis
from copy import deepcopy


def AddFitMetadata(path, h5_path, A, extra_dict=None):
        h5 = h5py.File(path)
        
        try: 
            h5.create_group(h5_path)
        except:
            pass

        fa = h5[h5_path].attrs
        fit = A.SaveFit()

        for key, val in fit.items():
            if key not in ['Data', 'energies', 'loglike', 'PSC']:
                fa.create('flux_'+key,val['flux'])
                fa.create('fluxunc_'+key,val['fluxunc'])        

        fa.create('loglike_total',np.sum(A.loglike))
        fa.create('loglike',A.loglike)
        fa.create('energies',A.central_energies)
        fa.create('bins', A.bin_edges)
        fa.create('irf', A.irf)
        fa.create('evclass', A.evclass)
        fa.create('convtype', A.convtype)
        fa.create('phfile', A.phfile)
        fa.create('tag', A.tag)

        if extra_dict is not None:
            for key, val in extra_dict.items():
                if key == 'residual':
                    try:
                        del h5[h5_path+'/residual']
                    except: 
                        pass

                    h5.create_dataset(h5_path+'/residual', data=val, dtype='float32')
                    print 'Saving new residual... Shape = ', val.shape

              #       try:
              #           del h5[h5_path+'/residual']
        		    # except:
        			   #  #pass
              #           h5.create_dataset(h5_path+'/residual', data=val, dtype='float32')
                else:
                    fa.create(key, val)
        h5.close()


def LoadModel(basedir, galprop_tag):
    # Load various diffuse models and run fits.
    print 'Running Analysis for model', galprop_tag
    A = Analysis.Analysis(tag='P7REP_CLEAN_V15_calore', basepath='/pfs/carlson/GCE_sys/')
    A.GenSquareMask(l_range=[-20.,20.], b_range=[-20.,20.], plane_mask=2.)
    A.BinPhotons(infile='binned_photons_'+A.tag+'.npy')
    # Load 2FGL 
    A.AddPointSourceTemplate(fixNorm=True,pscmap='PSC_3FGL_with_ext.npy')
    A.CalculatePixelWeights(diffuse_model='fermi_diffuse_'+A.tag+'.npy',psc_model='PSC_P7REP_CLEAN_V15_calore_fgl2.npy',
                        alpha_psc=5., f_psc=0.1)
    A.AddIsotropicTemplate(fixNorm=False, fixSpectrum=False) # External chi^2 used to fix normalization within uncertainties
    
    A.AddFermiBubbleTemplate(template_file='./bubble_templates_diskcut30.0.fits', 
                         spec_file='./reduced_bubble_spec_apj_793_64.dat', fixSpectrum=False, fixNorm=False)
    
    
    A.AddHDF5Template(hdf5file=basedir +'/'+ galprop_tag+'.hdf5',verbosity=1, multiplier=2., bremsfrac=1.25, 
                  E_subsample=2, fixSpectrum=False, separate_ics=False)
    return A



def Analyze(basedir, galprop_tag, A, analysis=0):

    if analysis == 0:
        #--------------------------------------------
        # GC fit without DM
        A.RunLikelihood(print_level=0, tol=2e2, precision=None, minos=True)[0]
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/GC_no_dm/', A=A, extra_dict={'residual':np.array([A.residual[i]*A.mask for i in range(A.n_bins)])})
        #AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/GC_no_dm/', A=A)
        

        #--------------------------------------------
        # GCE Fit
        A.ResetFit()    
        A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                       r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None,)
        A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=True)[0]

        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/GC_P7/', A=A, extra_dict={'residual':np.array([A.residual[i]*A.mask for i in range(A.n_bins)])})
        #AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/GC/', A=A)

    elif analysis == 1:
        #--------------------------------------------
        # Scan Slope
        gammas = np.linspace(.75,1.5,31)
        loglike_total, loglike, dm_spec, dm_spec_unc = [], [], [], []

        for i_g, gamma in enumerate(gammas):
            A.ResetFit()    
            print 'axes offset fitting completed:', i_g/float(len(gammas))
            A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=gamma, 
                           r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None,)
            A.RunLikelihood(print_level=0, tol=2e2, precision=None, minos=False)[0]
            loglike.append(A.loglike)
            loglike_total.append(np.sum(A.loglike))
            E, spec, specUnc = A.GetSpectrum('DM')
            dm_spec.append(spec)
            dm_spec_unc.append(specUnc)
        
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/scan_gamma_P7/', A=A, 
                       extra_dict={'gamma': gammas,
                                   'loglike':loglike,
                                   'loglike_total':loglike_total,
                                   'dm_spec':dm_spec,
                                   'dm_spec_unc':dm_spec})
    elif analysis == 2:
        #--------------------------------------------
        # Scan axes ratio
        ars = np.linspace(.6,2,21)
        loglike_total, loglike, dm_spec, dm_spec_unc = [], [], [], []
        for i_ar, ar in enumerate(ars):
            print 'axes offset fitting completed:', i_ar/float(len(ars))
            A.ResetFit()    
            A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                           r_s=20.0, axesratio=ar, offset=(0, 0), spec_file=None,)
            A.RunLikelihood(print_level=0, tol=2e2, precision=None, minos=False)[0]
            loglike.append(A.loglike)
            loglike_total.append(np.sum(A.loglike))
            E, spec, specUnc = A.GetSpectrum('DM')
            dm_spec.append(spec)
            dm_spec_unc.append(specUnc)
        
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/scan_axesratio_P7/', A=A, 
                       extra_dict={'axesratio': ars,
                                   'loglike':loglike,
                                   'loglike_total':loglike_total,
                                   'dm_spec':dm_spec,
                                   'dm_spec_unc':dm_spec},)

    elif analysis == 3:
        # #--------------------------------------------
        # # Scan longitude offset
        #lons = np.linspace(-90,90,61)
        lons = np.linspace(-45,45,31)
        loglike_total, loglike, dm_spec, dm_spec_unc, TS = [], [], [], [], []
        
        for i_l, lon in enumerate(lons):
            print 'lon offset fitting completed:', i_l/float(len(lons)), 'lon: ', lon

            A.ResetFit()
            A.templateList['Bubbles'].fixSpectrum = True
            A.templateList['Bubbles'].fixNorm = True
            A.GenSquareMask(l_range=[-20.+lon,20.+lon], b_range=[-20.,20.], plane_mask=2.)

            A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=False)[0]
            ll_nodm = np.sum(A.loglike)

            A.ResetFit()
            A.AddDMTemplate(profile='NFW', limits=[-100,100], decay=False, gamma=1.25, r_s=20.0, axesratio=1, offset=(lon, 0), spec_file=None,)

            A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=False)[0]

            loglike.append(A.loglike)
            TS.append(2*(ll_nodm-np.sum(A.loglike)))

            loglike_total.append(np.sum(A.loglike))
            E, spec, specUnc = A.GetSpectrum('DM')
            dm_spec.append(spec)
            dm_spec_unc.append(specUnc)

        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/scan_longitude_P7/', A=A, 
                       extra_dict={'longitudes': lons,
                                   'loglike':loglike,
                                   'loglike_total':loglike_total,
                                   'dm_spec':dm_spec,
                                   'dm_spec_unc':dm_spec,
                                   'TS': TS},)

    #--------------------------------------------
    # localize
    elif analysis == 4:
        lons = np.linspace(-1,1,21)
        fval = np.zeros((len(lons), len(lons)))
        
        for i_l, lon in enumerate(lons):
            for i_b, lat in enumerate(lons):
                print 'lat/lon fitting completed:', (len(lons)*i_l + i_b)/float(len(lons)**2)
                A.ResetFit()    
                A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                               r_s=20.0, axesratio=1, offset=(lon, lat), spec_file=None,)
                A.RunLikelihood(print_level=0, tol=2e2, precision=None, minos=False)[0]
                
                fval[i_b, i_l] = np.sum(A.loglike)
        
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/localize_P7/', A=A, 
                       extra_dict={'longitudes': lons,
                                   'latitudes': lons,
                                   'fval':fval},)

    
    elif analysis == 5:
        #--------------------------------------------
        # Scan Slope
        radius = np.linspace(2,28,14)
        loglike_total, loglike, dm_spec, dm_spec_unc = [], [], [], []
        A.AddIsotropicTemplate(fixNorm=False, fixSpectrum=False) # External chi^2 used to fix normalization within uncertainties

        
        A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                       r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None,)

        for i_r, r in enumerate(radius[:-1]):
            r1, r2 = r, radius[i_r+1]
                        
            mask = A.GenRadialMask(r1, r2, plane_mask=2, merge=False)
            # Now take NFW template and copy it, multiplied by the mask            
            A.templateList['ring_%i'%i_r] = deepcopy(A.templateList['DM'])    
            # Loop over energy and multiply by mask.
            for i_E in range(A.n_bins):
                A.templateList['ring_%i'%i_r].healpixCube[i_E] *= mask
        # Restore back to normal mask
        A.GenSquareMask(l_range=[-20.,20.], b_range=[-20.,20.], plane_mask=2.)
        A.DeleteTemplate('DM')
        A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=False)[0]
        
        r_bins = [(radius[i], radius[i+1]) for i in range(len(radius)-1)]
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/scan_radius_P7/', A=A, extra_dict={'radius':r_bins,})



    elif analysis == 7:
        radius = np.linspace(2,26,9)
        #radius = [2,4,10,15,20]
        loglike_total, loglike, dm_spec, dm_spec_unc = [], [], [], []
        A.AddIsotropicTemplate(fixNorm=False, fixSpectrum=False) # External chi^2 used to fix normalization within uncertainties
        A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                                   r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None)
        for i_r, r in enumerate(radius[:-1]):
            quadrants = [['E', -45, 45], ['N', 45, 135], ['W', 135, 225], ['S', 225, 315]]

            for quad, start_angle, stop_angle in quadrants:
                r1, r2 = r, radius[i_r+1]

                mask = A.GenRadialMask(r1, r2, plane_mask=2, merge=False, start_angle=start_angle, stop_angle=stop_angle)
                # Now take isotropic template and copy it, multiplied by the mask            
                A.templateList['ring_%i'%i_r + '_' + quad] = deepcopy(A.templateList['DM'])
                # Loop over energy and multiply by mask.
                for i_E in range(A.n_bins):
                    A.templateList['ring_%i'%i_r+ '_' + quad].healpixCube[i_E] *= mask
        # Restore back to normal mask
        A.GenSquareMask(l_range=[-20.,20.], b_range=[-20.,20.], plane_mask=2.)
        A.DeleteTemplate('DM')
        A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=False)[0]
        
        r_bins = [(radius[i], radius[i+1]) for i in range(len(radius)-1)]
        AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/quadrants/', A=A, extra_dict={'radius':r_bins,})

def AnalyzeMixture(basedir, galprop_tag, A, analysis, Pi0_Brems_PEB, Pi0_Brems_GAL):
    if analysis == 6:
        GalpropFracs = np.linspace(0.,1.,11)
        for i_GalpropFrac, GalpropFrac in enumerate(GalpropFracs):

            A.ResetFit() 
            A.templateList['Pi0_Brems'].healpixCube = Pi0_Brems_PEB*(1-GalpropFrac) + Pi0_Brems_GAL*GalpropFrac
            #--------------------------------------------
            # GC fit without DM
            A.DeleteTemplate('DM')
            A.RunLikelihood(print_level=0, tol=2e2, precision=None, minos=True)[0]
            AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/galprop_frac_no_dm_%1.2f/'%GalpropFrac, A=A, extra_dict=None)

            #--------------------------------------------
            # GCE Fit
            A.ResetFit()    
            A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.25, 
                           r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None,)
            A.RunLikelihood(print_level=1, tol=2e2, precision=None, minos=True)[0]
            AddFitMetadata(basedir +'/'+ galprop_tag+'.hdf5', h5_path='/fit_results/galprop_frac_%1.2f/'%GalpropFrac, A=A, extra_dict=None) 

import sys

if __name__ == "__main__":
    basedir, galprop_tag, analysis = sys.argv[1:4]
    if int(analysis) == 6:
        num = galprop_tag.split('_')[-2]
        print "working on", num, str(int(num)+28)
        A2 = LoadModel(basedir,galprop_tag.replace(num, str(int(num)+28)))
        Pi0_Brems_GAL = A2.templateList['Pi0_Brems'].healpixCube.copy()
        del(A2)
        A1 = LoadModel(basedir,galprop_tag)
        Pi0_Brems_PEB = A1.templateList['Pi0_Brems'].healpixCube.copy()
            
        AnalyzeMixture(basedir,galprop_tag, A1, int(analysis), Pi0_Brems_PEB, Pi0_Brems_GAL)    
    
    A = LoadModel(basedir,galprop_tag)
    Analyze(basedir,galprop_tag, A, int(analysis))





    

    #A.ResetFit()




# Run Analysis at GC
# Run Analysis without DM template. 
# Scan NFW slope
# Scan axis ratio
# scan offset. 
# Localize? 

