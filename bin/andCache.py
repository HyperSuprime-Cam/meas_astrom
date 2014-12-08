#!/usr/bin/env python

from lsst.meas.astrom.multiindex import AstrometryNetCatalog

def run(outName):
    AstrometryNetCatalog.fromEnvironment(False).writeCache(outName)
    print "Wrote cache %s" % outName

if __name__ == "__main__":
    import os

    outName = None
    if "ASTROMETRY_NET_DATA_DIR" in os.environ:
        outName = os.path.join(os.environ["ASTROMETRY_NET_DATA_DIR"], "andCache.fits")

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--outName", default=outName, required=(outName is None),
                        help="Output filename (default=$ASTROMETRY_NET_DATA_DIR/andCache.fits)")
    args = parser.parse_args()
    run(args.outName)
