# RFMIP-RRTMG

Piggy backs off of Robert Pincus's master branch of [RTE+RRTMGP](https://github.com/RobertPincus/rte-rrtmgp). Drivers are developed for RRTMG instead of RRTMGP, though, so that the timing of both can be compared to each other.

1. set RRTMGP_ROOT for scripts and Makefiles ()

```
export RRTMGP_ROOT=/rd47/scratch/RRTMGP/RRTMGP_Git
```

2. set TIMING_DIR, NCHOME, NFHOME for Makefiles (and build the libraries if they do not yet exist)

```
export NCHOME=/nas/project/p1770/dependencies
export NFHOME=/nas/project/p1770/dependencies
export TIME_DIR=/rd47/scratch/RRTMGP/RRTMGP_Git/examples/rfmip-clear-sky/GPTL-v5.5.3
```

3. clone branch *recursively*

```
git clone --recursive git@github.com:pernak18/rte-rrtmgp.git RRTMG
git checkout rrtmg
cd RRTMG
export RRTMG=`pwd`
```

4. build the RRTMG drivers

```
% pwd
${RRTMG}/rrtmg_LW
% make rrtmg_lw
% pwd
${RRTMG}/rrtmg_SW
% make rrtmg_sw
```

The build will generate `rrtmg_rfmip_lw` and `rrtmg_rfmip_sw` in their respective directories.

5. run the timing wrapper

```
RFMIP_wrapper.py -e $RRTMGP_ROOT/examples/rfmip-clear-sky/rrtmgp_SW/rrtmgp_rfmip_sw
RFMIP_wrapper.py -e $RRTMGP_ROOT/examples/rfmip-clear-sky/rrtmgp_LW/rrtmgp_rfmip_lw
RFMIP_wrapper.py -e $RRTMG/rrtmg_SW/rrtmg_rfmip_sw
RFMIP_wrapper.py -e $RRTMG/rrtmg_LW/rrtmg_rfmip_lw
```

Each run will generate 40 text files in the working directory (`TIMING.iter.block_size.model`, where for now iter is 1-10, block_size is [4, 8, 100, 1800], and model is either RRTMGP or RRTMG). For the purposes of step 6, let's assume that timing results have been organized into `rrtmg_lw`, `rrtmg_sw`, `rrtmgp_lw`, and `rrtmgp_sw` subdirectories.

6. Run the statistics script

```
% pwd
${RRTMG}/examples/rfmip-clear-sky
% stats_RRTMG_vs_RRTMGP.py -gl ${RRTMG}/rrtmg_LW -gs ${RRTMG}/rrtmg_SW -gpl ${RRTMG}/examples/rfmip-clear-sky/rrtmgp_LW -gps ${RRTMG}/examples/rfmip-clear-sky/rrtmgp_SW
```

The input directories for each model (RRTMG and RRTMGP) and domain (LW and SW) can be any name as long as the `timing.*.*.*` files from step 5 are in them.
