def rstparser(fh):
    for l in fh:
        if l.startswith("Bayes Empirical Bayes (BEB) probabilities for 4 classes (class)"):
            next(fh)
            next(fh)
            for l in fh:
                l = l.rstrip()
                if l == '':
                    return
                else:
                    f = l.split()
                    yield [ int(f[0]), f[1], float(f[2]), float(f[3]), float(f[4]), float(f[5]), int(f[7][0]) ]
