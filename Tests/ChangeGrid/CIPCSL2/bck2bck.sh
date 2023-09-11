#!/usr/bin/tcsh
foreach item (*.bck2)
   echo "I like $item"
end
foreach i (*.bck2)
mv $i $i:r.bck
end
