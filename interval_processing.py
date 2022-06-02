#!/usr/bin/env python3

"""
This is a python program that implements the task
"""

import sys, argparse

def usage():
    print("""Usage: python3 --input1 file1 --input2 file2 --length genomelength""")
    sys.exit(1)


def meanval(sample):
    try:
        result = round(sum(sample)/len(sample),4)
    except ZeroDivisionError as ze:
        result = str(e)
    except TypeError as te:
        result = str(te)
    return(result)
    

def pearsoncoef(sample1, sample2):
    if len(sample1) != len(sample2):
        result = str("Samples size not equal")
    elif len(sample1) == len(sample2):
        try:
            N = len(sample1)
            data = zip(sample1, sample2)
            numer = sum( (x-meanval(sample1))*(y-meanval(sample2)) for (x,y) in data)
            denom = (( sum((x-meanval(sample1)) ** 2 for x in sample1) ** 0.5) * ( sum((y-meanval(sample2)) ** 2 for y in sample2)** 0.5))
            result = round(numer/denom, 4)
        except ZeroDivisionError as ze:
            result = str(e)
        except TypeError as te:
            result = str(te)
    return(result)
    

def intervalIndex(intervalfile):
    with open(intervalfile) as file:
        intervalIndex = []
        for i, line in enumerate(file):
            start, end = line.split("\t")
            intervalRange = range(int(start),int(end.strip()))
            for p in intervalRange:
                intervalIndex.append(p)
    return(intervalIndex)


def intersect(list1, list2):
    set1 = set(list1)
    inter = list(set1.intersection(list2))
    inter.sort()
    return(inter) 


def createGenomeIndex(genomelength):
    genomePos = range(0, genomelength)
    genomeIndex = []
    for i in genomePos:
        genomeIndex.append(i)
    return(genomeIndex)


def parseMeasurements(file):
    with open(file) as file:
        measurements = []
        for line in enumerate(file):
            mes = float(line[1])
            measurements.append(mes)
    return(measurements)


def main():
    if len(sys.argv) < 3:
        usage()
    parser = argparse.ArgumentParser(description = 'process files')
    parser.add_argument('--input1', help = 'paths to input file1')
    parser.add_argument('--input2', help = 'paths to input file2')
    parser.add_argument('--length', type = int, help = 'length ot target genome')
    args = parser.parse_args()

    file1 = args.input1
    file2 = args.input2
    length = args.length

    # check input files type
    if file1.lower().endswith('.s') and file2.lower().endswith('.s'):
        try:
            # construct the genome index
            genome = createGenomeIndex(length)
            interval1 = intervalIndex(file1)
            interval2 = intervalIndex(file2)
            overlap = len(intersect(interval1, interval2))
            result = "the overlap of the input interval files is: " + str(overlap)
        except TypeError as te:
            result = str(te)        
    elif (file1.lower().endswith('.f') and file2.lower().endswith('.f')):
            try: 
                sample1 = parseMeasurements(file1)
                sample2 = parseMeasurements(file2)
                if len(sample2) != length or len(sample2) != length:
                    result = "incomplete measurements file"
                else: 
                    correlation = pearsoncoef(sample1, sample2)
                    result = "the Pearson correlation coefficient of the two samples is: " + str(correlation)
            except TypeError as te:
                result = str(te)   
    elif (file1.lower().endswith('.s') and file2.lower().endswith('.f')) or (file1.lower().endswith('.f') and file2.lower().endswith('.s')):
            try:
                if file1.lower().endswith('.s') and file2.lower().endswith('.f'):
                    interval = intervalIndex(file1)
                    sample = parseMeasurements(file2)
                elif file1.lower().endswith('.f') and file2.lower().endswith('.s'):
                    sample = parseMeasurements(file1)
                    interval = intervalIndex(file2)
                if len(sample) != length:
                    result = "incomplete measurements file"
                else:
                    part = []
                    for i in interval:
                        part.append(sample[int(i-1)])
                    average = meanval(part)
                    result = "the mean value of the measurements is: " + str(round(average, 4))
            except TypeError as te:
                result = str(te) 
    else:   result = "Wrong type of input files" 
    print(result)

if __name__ == "__main__":
    main()