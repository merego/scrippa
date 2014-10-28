#!/bin/bash

# I will put here all the images with correct size, density ecc.

echo "Convert images";

for i in *.ps;
do
name=$(basename $i .ps);
echo "Conversion of $name from ps to png 300 dpi, 3.25 in";
convert -flatten -resize 975 -density 118 -rotate 90 $name.ps $name.png;
echo " ";
done

