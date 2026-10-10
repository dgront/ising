SIZE="128"
ISING=$HOME/src.git/monte_carlo_sims/ising/target/release/ising 

mkdir -p frames

for t in 1.2 1.5 1.75 2.0 2.15 2.269 2.4 2.7 3.0 3.5; 
do 
    $ISING -l $SIZE --outer-cycles 10000 -t $t \
	    --outer-per-frame 10 --png-frames "frames/f_"$SIZE"_"$t"_" >en_$SIZE"_"$t;
done

