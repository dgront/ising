SIZE="1024"

cd frames
for t in 1.2 1.5 1.75 2.0 2.15 2.269 2.4 2.7 3.0 3.5
do

  i=0
  for f in f_$SIZE"_"$t"_"*.png
  do
    ln -s $f fr_$SIZE"_"$t"_"$(printf '%06d' $i).png
    i=$((i+1))
  done
  
  ffmpeg -framerate 30 \
      -i fr_$SIZE"_"$t"_"%06d.png \
      -c:v libx264 \
      -pix_fmt yuv420p \
      ising_$SIZE"_"$t".mp4"
done
cd ../

#f_512_2.15_9930.png
