# ps26_TestBeam

## Repository for Hall D Pair Spectrometer June 2026 TRD Test Beam

################################  
### Single Run Analysis Workflow

Workflow on JLab Gluon compute nodes:   
Log in to Gluon nodes with JLab computing account (must have 2FA)
```
ssh -XY [$USERNAME]@scilogin.jlab.org  
Password: [Pin+OTP]  
ssh hallgw  
Password: [Pin+OTP]  
ssh gluon[100-150]  
Password: [JLab CUE]  
```
After cloning repo, make directory links   
```
./make_gluon_links.sh  
```
And use executable files for run analysis   
```
./trdclass_ps26.sh [$RUNNUMBER] [$MAXNUMBEROFEVENTS] [$FIRSTEVENT]  
```
With this workflow, raw .evio data files have already been processed into .root files that now live in the `ROOT/` directory. These .root data files are analyzed with the `trdclass_ps26.C` analysis macro. Output from there is saved in the `RootOutput/ps26` directory. The output from this in the form of a .root TTree file, a .root file containing all monitoring histograms, and a PDF containing all of those histograms.  

###############################  
