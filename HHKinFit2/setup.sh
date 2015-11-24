#!/bin/bash

if [[ -z "$KINFIT2PATH" ]];
then
    export KINFIT2PATH=$PWD
    export LD_LIBRARY_PATH=$KINFIT2PATH:$LD_LIBRARY_PATH
    echo "Set KINFIT2PATH and updated LD_LIBRARY_PATH"

    if [[ $HOST == *"naf"*  ]];
    then
        module load git
        module load root
        module load gcc
        echo "Initialized git and ROOT via module (NAF usage)"
    else
        echo "Git and ROOT not explicitly initialized (non NAF usage)"
    fi
    
else
    echo "setup already done"
fi

