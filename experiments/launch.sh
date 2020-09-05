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
dataset="../../data/SR_div.csv"
seeds="../../seeds.csv"
G=40

for pop_sz in 50 100 300 500
do
    output_1="SR_div_test_${pop_sz}.csv"
    parameters_1="$dataset -p ${pop_sz} -G $G -t $seeds -g ramped_hh -f 0.000000001 -k 3 -d 7 -n 1 -o ${output_1} -M 0.05 -c .9"
#    ./sreg $parameters_1
done

pop_sz=500
output_2="SR_div_test_pc85_pm15.csv"
parameters_2="$dataset -p ${pop_sz} -G $G -t $seeds -g ramped_hh -f 0.000000001 -k 3 -d 7 -n 1 -o ${output_2} -M 0.15 -c .85"
#./sreg $parameters_2

for k in 5 7 10
do
    output_3="SR_div_test_k${k}.csv"
    parameters_3="$dataset -p ${pop_sz} -G $G -t $seeds -g ramped_hh -f 0.000000001 -k ${k} -d 7 -n 1 -o ${output_3} -M 0.3 -c .6"
#./sreg $parameters_3
done

k=5
output_4="SR_div_test_Eletism.csv"
parameters_4="$dataset -p ${pop_sz} -G $G -t $seeds -g ramped_hh -f 0.000000001 -k $k -d 7 -n 1 -o ${output_4} -M 0.3 -c .6 -E"
#./sreg $parameters_4


output_5="SR_div_test_Noise.csv"
parameters_4="$dataset -p ${pop_sz} -G $G -t $seeds -g ramped_hh -f 0.000000001 -k 5 -d 7 -n 1 -o ${output_5} -M 0.3 -c .6"
./sreg $parameters_4
