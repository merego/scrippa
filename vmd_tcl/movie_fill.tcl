#!/usr/bin/tclsh
# make movie tcl

#mol load pdb fort.91


color Display Background black


pbc set {302.777 302.777 302.777}
pbc box

rotate y by 45
rotate x by 30
scale by 1.3
set prot 200

mol delrep 0 top
#for {set i 0} {$i < $prot} {incr i} {
set i 0
set j 0
while {$i < $prot} {
  set i [expr {$i + 2}]
  set j [expr {$j + 1}]
  mol delrep 0 top
  mol selection "beta 1 to $i"
  mol addrep top
  set filename snap.[format "%04d" $j]
  render snapshot movie/$filename.tga
#          render Tachyon movie/$filename.tac
#          exec /home/mereghpo/my_lib/vmd/tachyon_LINUX  movie/$filename.tac -res 800 800 -o movie/$filename.tga
#          exec convert movie/$filename.tga movie/$filename.png

  exec convert movie/$filename.tga movie/$filename.png
  exec rm movie/$filename.tga
}
