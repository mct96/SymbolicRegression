import pandas as pd
import sys
import re

filename = sys.argv[1]
outfile = sys.argv[2]


with open(filename, "rt", encoding="utf-8") as rf:
    with open(outfile, "wt", newline="\n", encoding="utf-8") as wf:
        for line in rf:
            line = line.strip()
            line = re.sub("\s+", ",", line)
            wf.write(line + "\n")

data = pd.read_csv(open(outfile, "rt", encoding="utf-8"), header=None)
data.to_csv(open(outfile, "wt", encoding="utf-8"), sep="\t", float_format="%.3f", index=False, header=None)
