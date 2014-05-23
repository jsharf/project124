import argparse
import math

def main():
    parser = argparse.ArgumentParser(description='Baseline method of detecting\
            CNVs')
    parser.add_argument('rfile', type=str, nargs=1, help='Reference Genome')
    parser.add_argument('donorreads', type=str, nargs=1, help='SRS reads')
    parser.add_argument('numErrs', type=int, nargs=1, help='read err threshold')
    args = parser.parse_args()
    reference_genome = []
    with open(args.rfile[0], 'r') as ref:
        header = ref.readline()
        for i in ref:
            reference_genome.append(i)
    genome = ''.join(reference_genome)
    gIndex = {}
    # split genome up into segments of 50/2
    numAllowedErrs = args.numErrs[0]
    segmentSize = (50/numAllowedErrs)
    for i in range(0, len(genome), 50):
        for s in range(numAllowedErrs):
            seg = genome[i+(s)*segmentSize:i+(s+1)*segmentSize]
            if seg not in gIndex:
                gIndex[seg] = []
            gIndex[seg].append(i+(s-1)*segmentSize)
    # print("Indexed Reference genome. Now Sequencing...")
    reads = []
    with open(args.donorreads[0], 'r') as readfile:
        header = readfile.readline()
        for i in readfile:
            reads.append(i)
    coverage = [0 for idx in range(0,len(genome), segmentSize)];
    for x in reads:
        line = x.split('\t')
        read = line[1]
        for i in range(numAllowedErrs):
            seg = read[(i)*segmentSize:(i+1)*segmentSize]
            if seg in gIndex:
                pos = gIndex[seg]
                # print("FOUND SEG {0} AT POS('S): {1}".format(seg, pos))
                #if len(pos) == 1:
                coverage[pos[0]/segmentSize] += 1
    # now that I've made a coverage graph, time to detect CNV
    nonzeroes = filter(lambda x : x>0, coverage)
    tot = sum(nonzeroes)
    avg = float(tot)/len(nonzeroes)
    #var = float(sum([(x-avg)**2 for x in nonzeroes]))/len(nonzeroes)
    #stddev = math.sqrt(var)
    #print("Avg: " + str(avg))
    #print("Standard Deviation: " + str(stddev))
    # threshold for values which are at least 1 stddev above average
    thresholdedCoverage = [x > 1.5*avg for x in coverage] 
    in_cnv = False
    beginning_cnv = 0
    currCNVLen = 0 
    cnvlist = []
    for i in range(len(thresholdedCoverage)):
        if (not in_cnv) and thresholdedCoverage[i]:
            currCNVLen = 0
            in_cnv = True
            beginning_cnv = i
        if in_cnv and thresholdedCoverage[i]:
            currCNVLen += 1
        if in_cnv and not thresholdedCoverage[i]:
            in_cnv = False
            cnvlist.append((i, currCNVLen*segmentSize))
    print '\n'.join([str(x[0]) + ", " + str(x[1]) for x in cnvlist]) 


if __name__ == '__main__':
    main()
