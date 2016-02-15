import os
inputs_dirs = ["inputs", "inputs_tiny", "inputs_tiny_tiny"]
for inputs_dir in inputs_dirs:
    with open(os.path.join(inputs_dir, "rhos.tsv"), "r") as f_rhos:
        with open(os.path.join(inputs_dir, "wwph.suffixes"), "w") as f_suffixes:
            f_suffixes.write("RootNode:\n")
            for row in f_rhos:
                var, val = row[:-1].split("\t")
                val=str(max(float(val), 1e-6))
                f_suffixes.write("    {}:\n".format(var))
                f_suffixes.write("        CostForRho: {}\n".format(val))
