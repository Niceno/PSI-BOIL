#!/bin/tcsh

if ( $1 =~ "" || ($1 !~ "tec" && $1 !~ "visit") ) then
  echo "Wrong syntax!"
  echo "./verify.sh tec    -> use tecplot for visualization"
  echo "./verify.sh visit  -> use visit for visualization"
  exit
endif

echo $1

if ( $1 =~ "tec" ) then
  echo "use tecplot"
else
  echo "use visit"
endif

