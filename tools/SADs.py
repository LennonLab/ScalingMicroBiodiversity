import os
import sys

mydir2 = os.path.expanduser("~/Desktop/")

OrC = 'closed'

def GetSADsFromBiom_labeled(path, dataset):

    minS = 2

    IN = path + '/' + dataset + '-SSADdata.txt'
    n = sum(1 for line in open(IN))

    SiteDict = {}

    #print 'Starting build of SiteDict'
    with open(IN) as f:

        for d in f:

            #print 'Reading in SSAD data. Lines left:', n
            n -= 1

            if d.strip():

                d = d.split()
                #species = d[0]
                sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if sample not in SiteDict:

                        SiteDict[sample] = [abundance]

                    else:
                        SiteDict[sample].append(abundance)

    #print 'Finished building SiteDict'


    OUT = open(path + '/' + dataset + '-SADs.txt','w+')

    SADs = SiteDict.values()
    SiteDict = {}
    n = len(SADs)

    for i, sad in enumerate(SADs):
        print n - i

        if len(sad) >= minS:
            print >> OUT, sad

    OUT.close()

    print 'Finished generating SAD file'
    return

GetSADsFromBiom_labeled(mydir2 +'data/micro/EMP'+OrC, 'EMP'+OrC)
