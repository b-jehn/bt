import csv

minimizer_lengths = range(15,25)
thresholds = [16, 32, 64]

with open("runtime.lp", "w", encoding="utf-8") as f:
    f.write("Minimize\n")
    f.write(" obj: t_R_L1 + p_K_L1 t_l_L1 + (L1 - p_K_L1) (t_R_L2 + p_K_L2 t_l_L2 + (L1 - p_K_L2) (t_R_L3 + p_K_L3 t_l_L3 + (L1-p_K_L3) t_H))\n")
    f.write("Subject To\n")
    
    for level in range(1, 4):
        f.write(f"c{level}: p_K_L{level} =")
        for m in minimizer_lengths:
            with open(f'stats_m{str(m)}.csv') as mf:
                f.write(f" sm_{str(m)}_L{level}")
                reader = csv.reader(mf, delimiter=',')
                interestingrows=[row for idx, row in enumerate(reader) if idx in thresholds]
                for row in interestingrows:
                    t = row[0]
                    f.write(f" st_{t}_L{level} {row[2]}")
        f.write("\n")
    
    for level in range(1, 4):
        f.write(f"c{level+3}: 1 =")
        for m in minimizer_lengths:
            f.write(f" sm_{str(m)}_L{level}")
        f.write("\n")
        
    for level in range(1, 4):
        f.write(f"c{level+6}: 1 =")
        for t in thresholds:
            f.write(f" st_{t}_L{level}")
        f.write("\n")


