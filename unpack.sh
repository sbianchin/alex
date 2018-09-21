#!/bin/bash

cd ../dat2root

./run.sh $1


cd -

root -e Tree2Histo_Launcher.C\($1\)


