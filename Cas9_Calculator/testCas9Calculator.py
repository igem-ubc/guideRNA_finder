from Cas9_Calculator import *

if __name__ == "__main__":

    guideSequence = 'TACGTACACAAGAGCTCTAG'   
    Cas9Calculator = clCas9Calculator(['/home/connor/Bioinformatics/iGEM/NC_000913.gb'])
    sgRNA1 = sgRNA(guideSequence, Cas9Calculator)
    sgRNA1.run()
    # sgRNA1.exportAsDill()
    sgRNA1.printTopTargets()
