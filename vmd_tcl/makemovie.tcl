	# Set ambientocclusion
	display ambientocclusion on

	# get the number of frames in the movie
	set num [molinfo top get numframes]


	# loop through the frames
	for {set i 0} {$i < $num} {incr i} {
		# go to the given frame
		animate goto $i
		# take the picture
		set filename snap.[format "%04d" $i].tac
		render Tachyon movie/$filename
exec /home/mereghpo/my_lib/vmd/tachyon_LINUX  movie/$filename -res 800 800 -o movie/$filename.tga
exec convert movie/$filename.tga movie/$filename.ppm
	}

