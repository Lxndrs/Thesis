#!/bin/bash
echo "Starting multi fit protocol"
echo "Starting loop over folders"
for i in 00 01 02 03 04 05 06 07 08 09;
do 
    cd Multi-Fit/samp$i
    for j in LV2 LV2b LV3 LV4 LV5 177-341 169-338 180-331 168-328;
    do
	python ../../read-region-file-one-shock-fig.py --proplyd $j --region LV-positions-2018-will-samp$i.reg --tfit 180
    done
    cd ../..
done
echo "Finished!!"
