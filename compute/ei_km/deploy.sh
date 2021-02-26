#!/bin/bash

proj_root='/home/michael/Documents/Schol/Project/KM_Neuron/compute'

lyap=$(lab_util -host lyapunov)
poin=$(lab_util -host poincare)

rsync -av -e ssh --exclude='temp/' --exclude='data/' --exclude='deploy.sh' $proj_root'/ei_km' $lyap:Hall/
rsync -av -e ssh --exclude='temp/' --exclude='data/' --exclude='deploy.sh' $proj_root'/ei_km' $poin:Hall/

rsync -av -L -e ssh $proj_root'/km_util' $lyap:python_lib/
rsync -av -L -e ssh $proj_root'/km_util' $poin:python_lib/

