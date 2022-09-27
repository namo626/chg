# Fix ATCF storm files
from datetime import datetime as dt

#fname = 'bal012001.dat'
fname = 'allison.dat'
outfile = fname + '.fixed'

# Find the maximum RMW to fill the missing values
def findMaxRMW():
    maxRMW = 0
    with open(fname, 'r') as f1:
        while True:
            line = f1.readline();
            if not line:
                break

            # data is not missing
            if len(line) > 98:
                rmw = int(line.split(',')[19])
                if maxRMW < rmw:
                    maxRMW = rmw

    return maxRMW

def run():
    #maxRMW = findMaxRMW()
    with open(fname, 'r') as f1, open(outfile, 'w') as out:

        line = f1.readline()
        date_prev = dt.strptime(line.split(",")[2].strip(),"%Y%m%d%H")
        total_time = int(date_prev.hour)
        line = line[:30] + ("%03d" % total_time) + line[33:]
        # if missing data, append it (put RMW as 0)
        if len(line) <= 98:
                #line = '%s%5d,%5d,%4d%30d,%4d,%11s\n' %(line.strip(), 1013, -999, 0,-99,-99,name)
            line= '%s%5d,%5d,%4d\n' %(line.strip(), 1013, -999, 0)

        out.write(line)

        while True:
            line = f1.readline();
            if not line:
                break

            # fix the increment time (assume ADCIRC starts at hour 00)
            date_now = dt.strptime(line.split(",")[2].strip(),"%Y%m%d%H")
            inc_hours = int((date_now - date_prev).seconds/3600)
            total_time += inc_hours
            date_prev = date_now

            line = line[:30] + ("%03d" % total_time) + line[33:]

            # if missing data, append it (put RMW as 0)
            if len(line) <= 98:
                #line = '%s%5d,%5d,%4d%30d,%4d,%11s\n' %(line.strip(), 1013, -999, 0,-99,-99,name)
                line = '%s%5d,%5d,%4d\n' %(line.strip(), 1013, -999, 0)

            out.write(line)


    print("Done.")
