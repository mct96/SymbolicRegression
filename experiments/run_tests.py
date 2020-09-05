import os
import sys
import subprocess

path = "../bin/src/sreg"
dataset = "../data/SR_div.csv"
generation_method = "-g ramped_hh"
threshold = "-f 0.01"
k = "-k 10"
depth = "-d 7"
n_vars = "-n 1"
mutation = "-M 0.05"
crossover = "-c 0.9"
generations = "-G 90"
test = "-t ../seeds.csv"
parameters = [path, dataset, generation_method, threshold, k, depth, n_vars,
              mutation, crossover, generations, test]
# test population
for pop in [50, 100, 300, 500]:
    result = f"-o SR_div_{pop}.csv"
    cur_parameters = [path, dataset + " " + result  ]
    print(subprocess.run(cur_parameters, capture_output=True))
    
