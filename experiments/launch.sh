bin_dir=../bin

# create a build directory.
if [[ ! -e $bin_dir ]]; then
    makdir -p $bin_dir
else
    cd $bin_dir
fi

# start build.
cmake ..                        # generate build files.
cmake --build .                 # execute generators tools.
cd src                          # go to exetable generated path.

# start tests.
dataset="../../data/concrete.csv"
seeds="../../seeds.csv"

for pop_sz in 50 100 300 500
do
    output_1="test_pop_sz_$pop_sz_prob.csv"
    parameters_1="$dataset -p $pop_sz -G 40 -t $seeds -g ramped_hh -f 0.1 -k 3 -d 7 -n 8 -o $output_1 -M 0.3 -c .6 -E"
    ./sreg $parameters_1
done

#echo $parameters_1
#./sreg $parameters_1
