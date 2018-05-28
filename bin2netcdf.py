import os
import glob
from cdo import *
cdo = Cdo()

# Get a list of all files in the output directory
file_list = glob.glob('output/*')

for ifile in file_list:

    path, file = os.path.split(ifile)

    # Flux correction files
    if 'fluxcorr' in file:
        # Create ctl file
        ctl = open(path+'/'+'flx.sample.greb.ctl','w')
        ctl.write("""dset ^%s
        undef 9.e27
        xdef  96 linear 0 3.75
        ydef  48 linear -88.125 3.75
        zdef   1 linear 1 1
        tdef  730 linear 1jan0  12hr
        vars 3
        tsurf_corr  1 0 tsurf_corr
        tocean_corr  1 0 tocean_corr
        q_corr  1 0 q_corr
        endvars"""%(file))
        ctl.close()

        # Write netcdf
        cdo.import_binary(options='-f nc', input=path+'/'+'flx.sample.greb.ctl', output=path+'/'+file+'.nc')
        os.system('rm ' + path+'/'+'flx.sample.greb.ctl')

    # Control files
    if 'control' in file:
        # Create ctl file
        ctl = open(path+'/'+'control.sample.greb.ctl','w')
        ctl.write("""dset ^%s
        undef 9.e27
        xdef  96 linear 0 3.75
        ydef  48 linear -88.125 3.75
        zdef   1 linear 1 1
        tdef  60 linear 15jan0  1mo
        vars 9
        tsurf  1 0 tsurf
        tatm  1 0 tatm
        tocean  1 0 tocean
        q  1 0 q
        albedo 1 0 albedo
        precip 1 0 precip
        dqrain 1 0 dqrain
        dqeva 1 0 dqeva
        dqcrcl 1 0 dqcrcl
        endvars"""%(file))
        ctl.close()

        # Write netcdf
        cdo.import_binary(options='-f nc', input=path+'/'+'control.sample.greb.ctl', output=path+'/'+file+'.nc')
        os.system('rm ' + path+'/'+'control.sample.greb.ctl')

    # Scenario files
    if 'scenario' in file:
        # Create ctl file
        ctl = open(path+'/'+'scenario.sample.greb.ctl','w')
        ctl.write("""dset ^%s
        undef 9.e27
        xdef  96 linear 0 3.75
        ydef  48 linear -88.125 3.75
        zdef   1 linear 1 1
        tdef  60 linear 15jan0  1mo
        vars 9
        tsurf  1 0 tsurf
        tatm  1 0 tatm
        tocean  1 0 tocean
        q  1 0 q
        albedo 1 0 albedo
        precip 1 0 precip
        dqrain 1 0 dqrain
        dqeva 1 0 dqeva
        dqcrcl 1 0 dqcrcl
        endvars"""%(file))
        ctl.close()

        # Write netcdf
        cdo.import_binary(options='-f nc', input=path+'/'+'scenario.sample.greb.ctl', output=path+'/'+file+'.nc')
        os.system('rm ' + path+'/'+'scenario.sample.greb.ctl')
