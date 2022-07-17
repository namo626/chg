import numpy as np
import sys
import random
import argparse

parser = argparse.ArgumentParser(description="Create a copy of fort.13 with modified Manning's N")
parser.add_argument("fort13", help="Name of fort.13 file")
parser.add_argument("fort14", help="Name of fort.14 file")
parser.add_argument("--halve", help="Halve all Manning's N values", action="store_true")
parser.add_argument("--whole", help="Set domain of change to global", action="store_true")
args = parser.parse_args()

"""rand_manning.py

USAGE:
    $ python rand_manning.py <fort13> <fort14>

where fort13 and fort14 are the names of fort.13 and fort.14 files.
The output is written to a new file fort13.rand.
Requires numpy >= 1.16.0.

"""
#############################################################
# Main program

# Box to randomize manning's n
xmin = -94.
xmax = -93.
ymin = 29.
ymax = 30.

# Galveston bay
xmin = -96.
xmax = -94.25
ymin = 29.
ymax = 30.

def main():
    #fort13 = sys.argv[1]
    #fort14 = sys.argv[2]
    fort13 = args.fort13
    fort14 = args.fort14

    arr14 = load14(fort14)
    write13(fort13, arr14)


#############################################################
# Functions

# Load fort.14 into memory to access x and y coordinates for each node
def load14(fname):
    # Get the number of nodes
    with open(fname) as f:
        l1 = f.readline()
        l1 = f.readline()
        l1 = l1.split()
        NP = int(l1[1])
        print("NP = %d" %NP)

    arr = np.loadtxt(fname, skiprows=2, usecols=(1,2), max_rows=NP)
    return arr

# check if a node lies inside the box
def check_node(node, arr):
    x = arr[node-1, 0]
    y = arr[node-1, 1]
    if (x <= xmax and x >= xmin) and (y <= ymax and y >= ymin):
        return True
    if args.whole:
        return True 
    return False

# Load and copy fort.13 with edited section
def write13(fort13, arr14):
    with open(fort13, 'r') as f1, open(fort13 + ".rand", "w") as f2:

        # Read until first "manning's n" keyword
        wflg = False
        while True:
            line = f1.readline()
            lines = line.split()
            f2.write(line)
            if lines[0] == "mannings_n_at_sea_floor":
                if not wflg:
                    wflg = True
                else:
                    break

        # Randomize the manning's n
        l = f1.readline()
        mannings_node = int(l.split()[0])
        f2.write(l)

        print("Modifying manning's n...")

        rand_count = 0
        for i in range(mannings_node):
            line = f1.readline()
            lines = line.split()

            # if node in box, randomize/modify the friction
            node = int(lines[0])
            node_man = float(lines[1])
            if check_node(node, arr14):
                if arg.halve:
                    new_man = node_man / 2.0
                else:
                    new_man = random.uniform(0.02, 0.2)

                f2.write("%d %.6f\n" % (node, new_man))
                rand_count += 1

            else:
                f2.write(line)

        print("Finished randomizing...")
        # Copy the rest of the file
        while True:
            line = f1.readline()
            if not line:
                break

            f2.write(line)

        # Print how many nodes were randomized
        print(rand_count)


# check number of edits
def count_diff():
    count = 0
    with open(fort13, 'r') as f1, open("fort.13.rand", "r") as f2:
        while True:
            l1 = f1.readline()
            l2 = f2.readline()
            if (not l1) or (not l2):
                break
            if not l1 == l2:
                count += 1


    return count


main()
