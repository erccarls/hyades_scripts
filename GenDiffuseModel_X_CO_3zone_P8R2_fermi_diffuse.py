import Analysis
import cPickle as pickle
import Tools
import multiprocessing as mp 
import pyfits
import numpy as np
import h5py
import sys


def RunFit(A, nrings=9, limit_inner=None, fix_xco=False):
    #-----------------------------------------------------
    # Templates are now added so we fit X_CO
    
    import GammaLikelihood as like
    
    fval, res = [], []
    for key, t in A.templateList.items():
        # local 
        if key not in ['FermiDiffuse', 'Isotropic','pi0HIHII', 'ics', 'Bubbles']:
            t.fixNorm = True
            t.fixSpectrum= True
            t.limits = [0.0,200.]
            t.value=1.
        else: 
            t.fixNorm = False
            t.limits=[0., None]

            
    print 'Running Local Ring Fit...'
    A.GenSquareMask(l_range=[-180.,180], b_range=[-85.,85.], plane_mask=8.)
    m, R = A.RunLikelihood( print_level=1, precision=None, tol=1e3)[:2]
    fval.append(m.fval)
    res.append(R)

    print 'isotropic value:', A.templateList['Isotropic'].value
    
    # if fix_xco is False:
    #     vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
    #     print "X_CO adjustment (this is not XCO value, it is multiplier for galdef values of xco):", vals

    #-------------------------------------------------------------------
    # Now we have fit the local X_CO (fixed).  Next we fit the outer galaxy
    A.templateList['Isotropic'].fixNorm = True
    A.templateList['Bubbles'].fixNorm = True

    print 'Running Outer Rings Fit...'
    A.GenSquareMask(l_range=[-180,-80], b_range=[-8.,8.], plane_mask=0)
    A.GenSquareMask(l_range=[80,180], b_range=[-8.,8.], plane_mask=0, merge=True)

    m, R = A.RunLikelihood( print_level=1, precision=None, tol=1e3)[:2]
    fval.append(m.fval)
    res.append(R)


    #-------------------------------------------------------------------
    # Now we fit the inner galaxy X_CO.

    print 'Running Inner Rings Fit...'
    A.GenSquareMask(l_range=[-80,80], b_range=[-8.,8.], plane_mask=0)
    m, R = A.RunLikelihood( print_level=1, precision=None, tol=5e3)[:2]
    fval.append(m.fval)
    res.append(R)

    # if fix_xco is False:
    #     vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
    #     print "X_CO adjustment (this is not XCO value, it is multiplier for galdef values of xco):", vals

    return m, fval, res





def WriteHDF5(fname, basedir, tag, m , nrings=9, fix_xco=False):
    """
    Build the diffuse model according to the best fit parameters and writes model+metadata to an HDF5 file.
    
    :param fname: Output filename for the HDF5 file
    :param basedir: Directory with galprop output files
    :param tag: galprop tag
    :param basedir: iminuit object. 
    :param nrings: number of galprop rings.
    """
    modf = h5py.File(fname, 'w')

    if m is not None:
        if fix_xco is False:
            X_CO = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
        else: 
            X_CO = np.ones(nrings)
    #modf = h5py.File(fname, 'w')
    #try:
    # Generate Groups
    # template_group = modf.create_group("templates")
    # fit_group = modf.create_group("fit_results")


    # Get data dimensions
    tmp = pyfits.open(basedir+'/bremss_healpix_54_'+tag+'.gz')
    energies = tmp[2].data.field(0)
    tShape = (len(energies), tmp[1].data.shape[0])
    print tShape
    del tmp # free memory

    pi0     = modf.create_dataset("/templates/pi0", tShape, dtype='float32',compression="gzip")
    pi0_0   = modf.create_dataset("/templates/pi0_0", tShape, dtype='float32',compression="gzip")
    brem   = modf.create_dataset("/templates/brem", tShape, dtype='float32',compression="gzip")
    brem_0 = modf.create_dataset("/templates/brem_0", tShape, dtype='float32',compression="gzip")
    ics_opt = modf.create_dataset("/templates/ics_opt", tShape, dtype='float32',compression="gzip")
    ics_fir = modf.create_dataset("/templates/ics_fir", tShape, dtype='float32',compression="gzip")
    ics_cmb = modf.create_dataset("/templates/ics_cmb", tShape, dtype='float32',compression="gzip")
    modf.create_dataset("/templates/energies", data=energies, dtype='float32',compression="gzip")
    # Now fill in the templates one by one.
    # Add fit metadata.
    # Add galdef metadata.


    #---------------------------------------------------------------
    # Create Diffuse Template from fitting results.
    def ReadFits(fname, length):
        d = pyfits.open(fname)[1].data
        return np.array([d.field(i) for i in range(length)])

    if m is not None:
        for i_ring in range(1,nrings+1):
            print "Adding HI/HII ring", i_ring

            pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            if fix_xco is False:
                pi0[...] += X_CO[i_ring-1]*ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            else:
                pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

            brem[...] += m.values['pi0HIHII']*1.25*ReadFits(basedir+'/bremss_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            brem[...] += m.values['pi0HIHII']*1.25*ReadFits(basedir+'/bremss_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            if fix_xco is False:
                brem[...] += 1.25*X_CO[i_ring-1]*ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            else: 
                brem[...] += m.values['pi0HIHII']*1.25*X_CO[i_ring-1]*ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            
            if i_ring == 1:
                pi0_0[...] += pi0
                brem_0[...] += brem

        print 'Adding ICS'
        ics_opt[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_1_healpix_54_'+tag+'.gz', len(energies)).clip(0)
        ics_fir[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_2_healpix_54_'+tag+'.gz', len(energies)).clip(0)
        ics_cmb[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_3_healpix_54_'+tag+'.gz', len(energies)).clip(0)

    else:
        for i_ring in range(1,nrings+1):
            print "Adding HI/HII/H2 ring", i_ring

            pi0[...] += ReadFits(basedir+'/pi0_decay_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            pi0[...] += ReadFits(basedir+'/pi0_decay_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            pi0[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

            brem[...] += ReadFits(basedir+'/bremss_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            brem[...] += ReadFits(basedir+'/bremss_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            brem[...] += ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

            if i_ring == 1:
                pi0_0[...] += pi0
                brem_0[...] += brem

        print 'Adding ICS'
        ics_opt[...] += ReadFits(basedir+'/ics_isotropic_comp_1_healpix_54_'+tag+'.gz', len(energies)).clip(0)
        ics_fir[...] += ReadFits(basedir+'/ics_isotropic_comp_2_healpix_54_'+tag+'.gz', len(energies)).clip(0)
        ics_cmb[...] += ReadFits(basedir+'/ics_isotropic_comp_3_healpix_54_'+tag+'.gz', len(energies)).clip(0)

#     except:
#         modf.close()
    
    try: 
        modf.close()
        print 'Closed HDF5 file.' 
    except: 
        print 'Failed to close HDF5 file.' 
        pass
    return


def AddMetadata(fname, basedir, tag, A, m, fval=None, fix_xco=False):
    # Parse the galprop file into a dict.
    
    print 'Adding fit metadata'
    
    try: h5.close() 
    except: pass
    h5 = h5py.File(fname)
    
    
    if m is not None and A is not None:
        try:
            fit_results = h5.create_group("/fit_results/global")
        except: 
            fit_results = h5['/fit_results/global']
            
        fa = fit_results.attrs
        fa.create('globalvalues', m.values.items())
        fa.create('globalvaluesUnc', m.errors.items())
        fa.create('globalfval', m.fval)

        if fval is not None:
            fa.create('localfval', fval[0])            
            fa.create('outerfval', fval[1])            
            fa.create('innerfval', fval[2])            

        fa.create('globale_bins', A.bin_edges)
        fa.create('globalirf', A.irf)
        fa.create('globalevclass', A.evclass)
        fa.create('globalconvtype', A.convtype)
        fa.create('globalphfile', A.phfile)
        fa.create('globaltag', A.tag)
        h5.create_dataset('/fit_results/globalmask', data=A.mask, dtype='float32')
    h5.close()
    
try:
    modf.close()
    print 'Closed HDF5 file'
except:pass



if __name__ == "__main__":
    if (len(sys.argv) not in [4,5,6]) : 
        raise("Incorrect number of args: <galprop output dir> <galprop tag> <galdef dir> [limit_inner] [fix_xco]")
    
    basedir, tag, path = sys.argv[1:4]
    

    fname = basedir+'/'+tag+'_XCO_P8.hdf5'
    #fname = basedir+'/'+tag+'_XCO_P8_MS04.hdf5'


    # Load the analysis
    A = Analysis.Analysis(tag='P8R2_CLEAN_V6_calore', fglpath='/pfs/carlson/gll_psc_v16.fit',  
        templateDir='/home/carlson/pfs/Extended_archive_v15/Templates', basepath='/pfs/carlson/GCE_sys/')
    # A.GenPointSourceTemplate(pscmap=(A.basepath + '/PSC_all_sky_3fgl.npy'))
    # A.BinPhotons(outfile='binned_photons_all_sky.npy')
    #A.GenSquareMask(l_range=[-180.,180], b_range=[-40.,40.], plane_mask=1.)
    A.BinPhotons(infile='binned_photons_P8R2_CLEAN_V6_calore.npy')
    # Load 2FGL 
    A.AddPointSourceTemplate(fixNorm=True, pscmap=('PSC_3FGL_with_ext.npy'))
    A.CalculatePixelWeights(diffuse_model='fermi_diffuse_'+A.tag+'.npy',psc_model='PSC_3FGL_with_ext.npy',
                            alpha_psc=5., f_psc=0.1)
    A.AddIsotropicTemplate(fixNorm=False, fixSpectrum=True) # External chi^2 used to fix normalization within uncertainties
    #A.PopulateROI([0,0],radius=360, fix_radius=360., include_point=False)


    # OPEN THE Extended PSC file and add it to the template list. 
    A.templateList['PSCExt'] = pickle.load(open('PSCExt.pickle', 'rb'))
    A.AddFermiBubbleTemplate(template_file='./bubble_templates_diskcut30.0.fits',spec_file='./reduced_bubble_spec_apj_793_64.dat', fixSpectrum=True, fixNorm=False)

    #A.AddDMTemplate(profile='NFW', limits=[None,None], decay=False, gamma=1.26, 
    #                r_s=20.0, axesratio=1, offset=(0, 0), spec_file=None,)
    
    

    #--------------------------------------------------------------------
    # Run the analysis
    #A = GenDiffuse(A, basedir=basedir, tag=tag, verbosity=1, nrings=9, fix_xco=fix_xco)


    A.AddFermiDiffuseModel(path, infile=None, outfile=None, multiplier=1., fixSpectrum=True)
    
    m, fval, res = RunFit(A, nrings=9, limit_inner=None, fix_xco=False)
    #WriteHDF5(fname=fname, basedir=basedir, tag=tag, m=m, nrings=9, fix_xco=fix_xco)
    AddMetadata(fname,basedir=None, tag=tag, A=A, m=m, fval=fval, fix_xco=False)
    #---------------------------------------------------------------------

    # WriteHDF5(fname=fname, basedir=basedir, tag=tag, m=None, nrings=9)
    # AddMetadata(fname,basedir=galdefdir, tag=tag, A=None, m=None)    
