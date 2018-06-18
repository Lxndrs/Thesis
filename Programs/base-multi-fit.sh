#!/bin/bash
echo "Starting multi fit protocol"
echo "Starting loop over folders"

for j in LV2 LV2b LV3 LV4 LV5 177-341 169-338 180-331 168-328;
do
    python read-region-file-one-shock-fig.py --proplyd $j --region LV-positions-2018-will.reg --tfit 180
    mv LV-bowshocks-xyfancy-positionswill-$j.save ../saves/LV-bowshocks-xyfancy-base-positionswill-$j.save
done
echo "Finished!!"