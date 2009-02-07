# -*- python -*-
#
# Setup our environment
#
import glob, os.path, re, os
import lsst.SConsUtils as scons

env = scons.makeEnv("meas_astrom",
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/meas/astrom/trunk/SConstruct $",
                    [["boost", "boost/version.hpp", "boost_system:C++"],
                     ["boost", "boost/version.hpp", "boost_filesystem:C++"],
                     ["boost", "boost/regex.hpp", "boost_regex:C++"],
                     ["python", "Python.h"],
                     ["m", "math.h", "m", "sqrt"],
                     ["cfitsio", "fitsio.h", "cfitsio", "ffopen"],
                     ["wcslib", "wcslib/wcs.h", "wcs"],
                     ["xpa", "xpa.h", "xpa", "XPAPuts"],
                     ["pex_exceptions", "lsst/pex/exceptions.h", "pex_exceptions:C++"],
                     ["utils", "lsst/utils/Utils.h", "utils:C++"],
                     ["daf_base", "lsst/daf/base.h", "daf_base:C++"],
                     ["pex_policy", "lsst/pex/policy/Policy.h", "pex_policy:C++"],
                     ["pex_logging", "lsst/pex/logging/Trace.h", "pex_logging:C++"],
                     ["security", "lsst/security/Security.h", "security:C++"],
                     ["daf_persistence", "lsst/daf/persistence/Persistence.h", "daf_persistence:C++"],
                     ["daf_data", "lsst/daf/data/LsstBase.h", "daf_data:C++"],
                     ["astrometry_net", "pnpoly.h", "pthread backend"], 
                     ["afw", "lsst/afw/image/MaskedImage.h", "afw"],
                    ]
                   )

env.libs["meas_astrom"] +=  env.getlibs("daf_base daf_data daf_persistence pex_logging pex_exceptions pex_policy security afw boost utils wcslib astrometry_net")
#
# Build/install things
#
for d in Split("doc examples include/lsst/meas/astrom/net python/lsst/meas/astrom/net lib src/net tests"):
    SConscript(os.path.join(d, "SConscript"))

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

if False:
    Alias("install", env.Install(env['prefix'], "python"))
    Alias("install", env.Install(env['prefix'], "include"))
    Alias("install", env.Install(env['prefix'], "lib"))
    #Alias("install", env.Install(env['prefix'], "pipeline"))
    Alias("install", env.Install(env['prefix'] + "/bin", glob.glob("bin/*.py")))
    Alias("install", env.InstallEups(env['prefix'] + "/ups", glob.glob("ups/*.table")))

    scons.CleanTree(r"*~ core *.so *.os *.o")
    #
    # Build TAGS files
    #
    files = scons.filesToTag()
    if files:
        env.Command("TAGS", files, "etags -o $TARGET $SOURCES")


env.Declare()
env.Help("""
LSST FrameWork packages
""")

