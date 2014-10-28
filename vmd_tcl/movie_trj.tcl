set i 0
set num [molinfo top get numframes]
while {$i < $num} {
  set i [expr {$i + 1}]
# go to the given frame
  animate goto $i
  set filename trj.[format "%04d" $i]
  render snapshot movie/$filename.tga
#          render Tachyon movie/$filename.tac
#          exec /home/mereghpo/my_lib/vmd/tachyon_LINUX  movie/$filename.tac -res 800 800 -o movie/$filename.tga
#          exec convert movie/$filename.tga movie/$filename.png

  exec convert movie/$filename.tga movie/$filename.png
  exec rm movie/$filename.tga
}
