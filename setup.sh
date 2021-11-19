#!/bin/sh

if [[ $1 == dg ]]; then
  cp ~/chg/dg-setup/* .
  cp ~/bin/adcprep-dg .
  cp ~/bin/dgswem .
  cp ~/bin/adcpost .
else
  cp ~/chg/adc-setup/* .
  cp ~/bin/adcprep .
  cp ~/bin/padcirc .
fi

cp ~/chg/Figuregen/Default2.pal .
cp ~/chg/Figuregen/shin.inp .
