import csv

minimizer_lengths = range(15,25)
thresholds = [16, 32, 64]

with open("runtime.lp", "w", encoding="utf-8") as f:
    #optimization goal:
    f.write("Minimize\n")
    f.write(" obj: t_R_L1 + p_K_L1 t_l_L1 + (1 - p_K_L1) (t_R_L2 + p_K_L2 t_l_L2 + (1 - p_K_L2) (t_R_L3 + p_K_L3 t_l_L3 + (1-p_K_L3) t_H))\n")
    
    #constraints:
    f.write("\nSubject To\n")
    
    for level in range(1, 4):
        f.write(f"c{level}: p_K_L{level} = ")
        for m in minimizer_lengths:
            with open(f'stats_m{str(m)}.csv') as mf:
                f.write(f"sm_{str(m)}_L{level} (")
                reader = csv.reader(mf, delimiter=',')
                for idx, row in enumerate(reader):
                    if idx not in thresholds:
                        continue
                    t = row[0]
                    f.write(f"{row[2]} st_{t}_L{level}")
                    if idx != thresholds[-1]:
                        f.write(f" + ")
                f.write(f")")
            if m != minimizer_lengths[-1]:
                f.write(f" + ")
        f.write("\n")
    
    for level in range(1, 4):
        f.write(f"c{level+3}: 1 = ")
        for m in minimizer_lengths:
            f.write(f"sm_{str(m)}_L{level}")
            if m != minimizer_lengths[-1]:
                f.write(f" + ")
        f.write("\n")
        
    for level in range(1, 4):
        f.write(f"c{level+6}: 1 = ")
        for t in thresholds:
            f.write(f"st_{t}_L{level}")
            if t != thresholds[-1]:
                f.write(f" + ")
        f.write("\n")
        
    #binary variables:
    f.write("\nBINARY\n")
    for level in range(1, 4):
        for m in minimizer_lengths:
            f.write(f"sm_{str(m)}_L{level} ")
        for t in thresholds:
            f.write(f"st_{str(t)}_L{level} ")
        f.write("\n")
        
    


