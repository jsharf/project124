from random import *
import argparse
import os
from time import time

rndgen = Random(time())

def get_random_genome(charset, length):
    global rndgen
    genome = []
    for i in range(length):
        genome.append(choice(charset))
    return genome

def mutate(genome, snps, pmutate, pcopy, lenCNV):
    global rndgen
    nextgen = list(genome)
    for i in snps:
        if (rndgen.uniform(0, 1) < pmutate):
            nextgen[i[0]] = i[2]
    segment_length = lenCNV
    segments = [nextgen[x:x+segment_length] for x in range(0, len(nextgen), segment_length)]
    addons = []
    for x in segments:
        if(rndgen.uniform(0, 1) < pcopy):
            addons.append(x)
    segments += addons
    nextgen = [inner for outer in segments for inner in outer]
    return nextgen 

def getReads(donor):
    rng = range(0, len(donor), 50)
    print len(donor)
    reads = [(seg, donor[seg:seg+50]) for seg in rng]
    shuffle(reads)
    return reads

def main():
    parser = argparse.ArgumentParser(description='Generate a reference and\
     donor genome')
    parser.add_argument('length', type=int, nargs=1, help='the length of the genome to generate')
    parser.add_argument('nSnps', type=int, nargs=1, help='the number of SNPs to create')
    parser.add_argument('nGens', type=int, nargs=1, help='the number of generations to simulate')
    parser.add_argument('pm', type=float, nargs=1, help='the probability of an snp mutation')
    parser.add_argument('pCopy', type=float, nargs=1, help='Probability of a segment copying in a generation')
    parser.add_argument('lenCNV', type=int, nargs=1, help='the size of CNV to simulate')
    args = parser.parse_args()
    global rndgen
    print "Genome length: " + str(args.length)
    print "Number SNPs: " + str(args.nSnps)
    print "Number Generations: " + str(args.nGens)
    print "Probability of Mutation: " + str(args.pm)
    print "Probability of Copy: " + str(args.pCopy)
    print "Length of CNV: " + str(args.lenCNV)
    charset = "ATGC"
    genome = get_random_genome(charset, args.length[0])
    snps = []
    for i in range(args.nSnps[0]):
        idx = rndgen.randint(0, len(genome))
        val = genome[idx]
        alts = (set(charset) - set(val))
        newbase = sample(alts, 1)[0]
        snp = (idx, val, newbase)
        snps.append(snp)
    donor = genome
    for i in range(args.nGens[0]):
        donor = mutate(donor, snps, args.pm[0], args.pCopy[0], args.lenCNV[0]);
    rfilestr = "ref"
    num = 0;
    rfilesuf = ".fasta"
    while (os.path.isfile(rfilestr+str(num) + rfilesuf)):
        num+=1
    glines = '\n'.join([''.join(genome[i:i+50]) for i in range(0, len(genome), 50)])
    with open(rfilestr+str(num)+rfilesuf, "w") as refFile:
        refFile.write(">Reference Genome\n"+glines)

    reads = getReads(donor)
    lines = [str(i) + "\t" + ''.join(read[1]) + "\t100\t" + str(read[0]) for i,read in enumerate(reads)]
    header = "id\tread_sequence\tinsert_length\tmapping_pos\n"
    dfilestr = "donor"
    dfilesuf = ".reads"
    num = 0
    while (os.path.isfile(dfilestr + str(num) + dfilesuf)):
        num+=1
    with open(dfilestr+str(num)+dfilesuf, "w") as donorFile: 
        donorFile.write(header + '\n'.join(lines));
    print "DONE!";

if __name__ == '__main__':
    main()
