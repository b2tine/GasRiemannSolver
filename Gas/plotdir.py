import matplotlib.pyplot as plt
import sys
import os


outdir = sys.argv[1].rstrip("/")

dirlist = list()
for name in os.listdir(outdir):
    dirname = os.path.join(outdir,name)
    if os.path.isdir(dirname):
        dirlist.append(dirname)

for subdir in dirlist:
    pngdir = os.path.join(subdir,"png")
    try:
        os.mkdir(pngdir)
    except FileExistsError:
        pass

    for filename in os.listdir(subdir):
        X = list()
        Y = list()
        
        if (filename.endswith(".txt")):
            for lines in open(os.path.join(subdir,filename),'r'):
                values = lines.split()
                X.append(float(values[0]))
                Y.append(float(values[1]))
            
            plt.plot(X,Y)
            plt.savefig(os.path.join(pngdir,filename.rstrip(".txt") + ".png"))
            plt.clf()
