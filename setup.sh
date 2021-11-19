#!/bin/sh

if [[ $1 == dg ]]; then
  cp ~/chg/dg-setup/* .
else
  cp ~/chg/adc-setup/* .
fi

cp ~/chg/Figuregen/Default2.pal .
cp ~/chg/Figuregen/shin.inp .
