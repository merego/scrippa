#!/usr/bin/tclsh
# Load clusters centroid generated using amber
# and create snapshot of the first N centroid  
# USAGE: vmd -dispdev text < make_cluster_snaps -args <parameterfile> <centroid-trajectory> <number-of-centroids-to-snap>
# EXAMPLE: vmd -dispdev text < vmd_select.tcl -args ../RNAFold.solv.prmtop clusters.mdcrd 3
#
# AUTHOR: Paolo Mereghetti, 2012
#
puts "Number of arguments $argc"
if { $argc != 4 } {
   puts "USAGE:  vmd -dispdev text < make_cluster_snaps -args <parameterfile> <centroid-trajectory> <number-of-centroids-to-snap>"
   puts "EXAMPLE: vmd -dispdev text < vmd_select.tcl -args ../RNAFold.solv.prmtop clusters.mdcrd 3"
   exit
}
set parameterfile [lindex $argv 0]
set trajectory [lindex $argv 1]
set num [lindex $argv 2]
#replace _ with spaces
#regsub -all {_} $vmdsel " " vmdsel
puts "###########################################"
puts "Parameter file :     $parameterfile"
puts "Trajectory:     $trajectory"
puts "Number of centroids to snap: $num"
puts "###########################################"

display resize 800 800
display projection orthographic
color Display Background white
axes location off

# Load amber parameter file
# and amber trajectory (ASCII format with perdiodic box)
#mol load parm7 ../RNAFold.solv.prmtop crdbox clusters.mdcrd
mol load parm7 $parameterfile crdbox $trajectory
display resetview

# Clear all
mol delrep 0 top
# Define representation
mol selection nucleic
mol representation newribbons
mol material AOEdgy
mol addrep top
# Define another representation
mol selection nucleic
mol representation hbonds 3.0 20.0 5.0
mol color colorID 16
mol addrep top


# Set ambientocclusion
display ambientocclusion on
display shadows on

# get the number of frames in the movie
#set num [molinfo top get numframes]


# loop through the frames
#set num 3
for {set i 0} {$i < $num} {incr i} {
# go to the given frame
animate goto $i
# Center and adjust zoom
display resetview
puts $i
scale by 2.0
# take the picture
set filename cluster.[format "%04d" $i]
render Tachyon cluster_snaps/$filename.tac
exec /usr/local/lib/vmd/tachyon_LINUXAMD64  cluster_snaps/$filename.tac -res 800 800 -o cluster_snaps/$filename.tga
exec convert cluster_snaps/$filename.tga cluster_snaps/$filename.png
#exec convert movie/$filename.tga movie/$filename.ppm
}

