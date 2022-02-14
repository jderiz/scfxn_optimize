
def parse(cpu_string):
    """
    This function parsese $SLURM_CPUS_PER_NODE into a comma seperated string of 
    cpus per node. e.g ...,18(x3),... --> ...,18,18,18,...
    """
    splits = cpu_string.split(",")
    res = ""

    for split in splits:
        split.split("(")

        if len(split) <= 1:
            res += split  # append single to result

            continue
        else:
            num_procs = split[0]
            repeats = split[1].strip(")").strip("x")

            for i in range(int(repeats)):
                res += num_procs+" "
    return res
