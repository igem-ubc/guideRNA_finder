#!/usr/bin/env python

import csv
import time
import datetime
import scipy.io
import math
import argparse
import re
import sys
from itertools import product
from time import time
from Bio import SeqIO

# import operator
# import dill


def get_options():
    """
    Function that uses argparse to parse the command-line arguments and return an object with the parameters
    :return: args - an object containing parameters parsed from the command-line
    """
    parser = argparse.ArgumentParser(description="A python script to look for off-target guide RNA binding sites using"
                                                 "a potential guide RNA sequence and a FASTA or Genbank (.gbk) file")
    parser.add_argument("-q", "--query",
                        required=True,
                        nargs="+",
                        help="File(s) containing sequences that need to be analyzed for off-target guideRNA binding.")
    
    parser.add_argument("-t", "--target_sequence",
                        required=True,
                        help="The gene sequence to be targetted")

    parser.add_argument("-m", "--model",
                        required=True,
                        help="The matlab model to use.")
    args = parser.parse_args()
    return args


def mers(length):
    """
    Generates multimers for sorting through list of k-mers based on user
    specification. Multimers generated act as the keys for generating a
    hashtable to eliminate undesired sequence patterns from those k-mers not
    found in the genome.

    Usage: mers(k) = 4^(k) unique k-mers
    """
    # Scales equally well as the old code, but simpler
    seq_list = list()
    nucleotides = ['A', 'T', 'C', 'G']
    all_kmers = product(nucleotides, repeat=length)
    for mer in all_kmers:
        seq_list.append(''.join(mer))
    return seq_list


def get_nggs(target_sequence):
    """

    :param target_sequence:
    :return:
    """
    target_seq = ""
    with open(target_sequence, 'r') as target:
        for line in target:
            line = line.upper()
            if line[0] == '>':
                pass
            else:
                target_seq += line.strip()

    nggs = list()
    length = 23
    counter = 0
    ngg_re = re.compile(r'.*GG$')
    while counter < len(target_seq) - length:
        candidate = target_seq[counter:counter+length]
        if re.match(ngg_re, candidate):
            nggs.append(candidate[:20])
        counter += 1

    sys.stdout.write("done.\n")
    return nggs


def identify_mer_positions(full_sequence, empty_mers, length=10):
    """
    Saves list of nucleotide positions in genome that all match a unique N-mer
    sequence. Counting begins at _ending_ of MER.

    Usage:   genomePositionsAtMers[mer_sequence] is a list of nucleotide positions
    within the inputted full_sequence that match mer_sequence
             If mer_sequence ends in a PAM site, then this can be used to match
             the first N-3 nt of a guide strand plus a PAM site sequence.
    """

    counter = 0
    while counter < (len(full_sequence)-length):
        word = full_sequence[counter:counter+length]
        if re.search('N', word):
            pass
        else:
            try:
                empty_mers[word].append(counter+length)
            except KeyError:
                sys.exit(word)
        counter += 1
    return empty_mers


def identify_target_sequence_matching_pam(pam_seq, positions_at_mers, full_sequence, target_sequence_length=20):
    """ Generates a list of target nucleotide sequences and corresponding nt
    positions for an inputted sequence that matched the pam_seq.
        Uses the positionsAtMers dictionary to accelerate the identification.
        Good for large genomes.

    Usage:  listOfTargets = identify_target_sequence_matching_pam('CGG', positionsAtMers, genome_sequence)
    """
    target_sequence_list = []
    all_mers = positions_at_mers.keys()
    mer_length = len(all_mers[0])
    list_of_mers_with_pam = [mer + pam_seq for mer in mers(mer_length - len(pam_seq))]
    for mer_with_PAM in list_of_mers_with_pam:
        nt_list = positions_at_mers[mer_with_PAM]
        for nt in nt_list:
            begin = nt-target_sequence_length - len(pam_seq)
            end = nt - len(pam_seq)
            if begin > 0 and end < len(full_sequence):  # Does not account for circular DNAs
                target_sequence = full_sequence[begin:end]
                target_sequence_list.append((target_sequence, nt))
    return target_sequence_list


class sgRNA(object):

    def __init__(self, guide_sequence, Cas9Calculator):

        self.guide_sequence = guide_sequence
        self.Cas9Calculator = Cas9Calculator

        self.partition_function = 1
        self.targetSequenceEnergetics = {}

        self.debug = False

    def run(self):

        begin_time = time()

        targetDictionary = self.Cas9Calculator.targetDictionary
        for (source, targets) in targetDictionary.items():
            self.targetSequenceEnergetics[source] = {}
            for fullPAM in self.Cas9Calculator.returnAllPAMs():
                dG_PAM = self.Cas9Calculator.calc_dG_PAM(fullPAM)
                dG_supercoiling = self.Cas9Calculator.calc_dG_supercoiling(sigmaInitial=-0.05, targetSequence=20 * "N")  #only cares about length of sequence
                for (targetSequence, targetPosition) in targetDictionary[source][fullPAM]:
                    dG_exchange = self.Cas9Calculator.calc_dG_exchange(self.guide_sequence, targetSequence)
                    dG_target = dG_PAM + dG_supercoiling + dG_exchange
                    self.targetSequenceEnergetics[source][targetPosition] = {'sequence': targetSequence, 
                                                                             'dG_PAM': dG_PAM, 
                                                                             'full_PAM': fullPAM, 
                                                                             'dG_exchange': dG_exchange,
                                                                             'dG_supercoiling': dG_supercoiling,
                                                                             'dG_target': dG_target}
                    self.partition_function += math.exp(-dG_target / self.Cas9Calculator.RT)

                    if self.debug:
                        print "targetSequence : ", targetSequence
                        print "fullPAM: ", fullPAM
                        print "dG_PAM: ", dG_PAM
                        print "dG_supercoiling: ", dG_supercoiling
                        print "dG_exchange: ", dG_exchange
                        print "dG_target: ", dG_target
                        print "Partition function (so far): ", self.partition_function

        end_time = time()

        print "Elapsed Time: ", end_time - begin_time
        print "Target Sequence: ", self.guide_sequence

    def printTopTargets(self, num_targets_returned=100):


        for (source, targets) in self.targetSequenceEnergetics.items():
            print "SOURCE: %s" % source


            # sort the targets items by the dG_target field
            sortedTargetList = sorted(targets.items(), key=lambda (k, v): v['dG_target'])  # sort by smallest to largest dG_target

            #print the index of the table
            print "POSITION\t\tTarget Sequence\t\tdG_Target\t\t% Partition Function"

            # for every element (0 -to- end of list) in the now sorted list of target sites found: do something
            for (position, info) in sortedTargetList[0:num_targets_returned]:

                #calculate the percent partition (Partision Function) for this target site
                percentPartitionFunction = 100 * math.exp(-info['dG_target'] / self.Cas9Calculator.RT) / self.partition_function

                #print the info to the screen: position, info["sequence"], info[dG_target] and the partition function we just calculated
                print "%s\t\t\t%s\t\t\t%s\t\t\t%s" % (str(position), info['sequence'], str(round(info['dG_target'], 2)), str(percentPartitionFunction) )

    def getResults(self):
        """
        This function should return the data that this run of the function calculated
        return an array of the best in the form[Target location, Partition Function] FOR THIS RUN OF THE CODE
        """
        # TODO: return the results of this run of the calculator
        # TODO: note: the above function does almost excactly what we want, we just need to pullout the first thing it prints out and return the first one of the sorted list

        output = [str("SOME INFO"), str("MORE INFO"), str("LAST BIT OF INFO")]

        for (source, targets) in self.targetSequenceEnergetics.items():
            print "SOURCE: %s" % source
            print "1"
            # sort the targets items by the dG_target field
            sortedTargetList = sorted(targets.items(), key=lambda (k, v): v['dG_target'])  # sort by smallest to largest dG_target

            top_hit = sortedTargetList[0]
            position = str(top_hit[0])  # assuming position is in index 0
            target_sequence = top_hit[1]['sequence']
            dg_target = str(round(top_hit[1]['dG_target'], 2))
            percentPartitionFunction = 100 * math.exp(
                -top_hit[1]['dG_target'] / self.Cas9Calculator.RT) / self.partition_function
            partition_function = str(percentPartitionFunction)
            output = [str(self.guide_sequence), position, target_sequence, dg_target, partition_function]

            return output


class clCas9Calculator(object):

    def __init__(self, filename_list, model_name='InvitroModel.mat', quickmode=True):

        self.quickmode = quickmode
        self.model_name = model_name
        data = scipy.io.loadmat(self.model_name)
        self.weights = data['w1']
        self.decNN = data['decNN']
        self.RT = 0.61597

        # the PAMs with the highest dG, ignoring other PAM sequences by setting their dG to 0
        self.PAM_energy = {'GGA': -9.8, 'GGT': -10, 'GGC': -10, 'GGG': -9.9, 'CGG': -8.1, 'TGG': -7.8, 'AGG': -8.1,
                           'AGC': -8.1, 'AGT': -8.1, 'AGA': -7.9, 'GCT': -7.1, 'GCG': -6.9, 'ATT': -7.2, 'ATC': -6.4,
                           'TTT': -7.6, 'TTG': -6.8, 'GTA': -7.4, 'GTT': -7.9, 'GTG': -7.7, 'AAT': -7, 'AAG': -7,
                           'TAT': -7.2, 'TAG': -7.2, 'GAA': -7.2, 'GAT': -7.3, 'GAC': -7.2, 'GAG': -7.3}

        self.init_target_finder(filename_list)

    def returnAllPAMs(self):

        for (PAMpart, energies) in sorted(self.PAM_energy.items(), key=lambda x: x[1]):  # PAMpart will be 'GGT'
            for nt in ('A', 'G', 'C', 'T'):  # nt + PAMpart will be all possible 'NGGT'
                yield nt + PAMpart

    def init_target_finder(self, filename_list, length=10):

        target_dictionary = {}

        # Create a list of all possible k-mers
        all_possible_mers = mers(length)
        # Search through the genome and add nucleotide positions for match to an N-mer
        empty_mers = {}
        for mer in all_possible_mers:
            empty_mers[mer] = []

        print "Number of k-mers: ", len(empty_mers.keys())

        full_sequence = ""
        for filename in filename_list:
            extension = filename.split('.')[-1]
            handle = open(filename, 'r')
            if extension in ['fa', 'fna', 'faa', 'fasta']:
                # File is in FASTA format
                fa_records = SeqIO.parse(handle, "fasta")
                for record in fa_records:
                    full_sequence += str(record.seq)
            elif extension in ['gb', 'gbff']:
                # File is a GenBank file
                genbank_record_obj = SeqIO.parse(handle, "genbank")
                record = genbank_record_obj.next()
                full_sequence += str(record.seq)
            else:
                sys.stderr.write("ERROR: Unable to detect file-type of " + filename + " by extension! Exiting...\n")
            handle.close()

        positions_at_mers = identify_mer_positions(full_sequence, empty_mers, length)
        query_files = ','.join(filename_list)
        target_dictionary[query_files] = dict()
        targetSequenceList = []
        for fullPAM in self.returnAllPAMs():
            targetSequenceList = identify_target_sequence_matching_pam(fullPAM, positions_at_mers, full_sequence)
            target_dictionary[query_files][fullPAM] = targetSequenceList
        self.targetDictionary = target_dictionary

        empty_mers.clear()

    def print_model_info(self):
        m = 0
        s = 0
        negative_val = 0
        for i, l in enumerate(self.decNN):
            for j, e in enumerate(l):
                if float(e) < 0:
                    negative_val += 1
                if i != j:
                    s += float(e)
                    m += 1
        meanNN=float(s)/float(m)

        sw = 0
        for w in self.weights:
            sw += w

        meanw = sw/len(self.weights)
        print 'average mismatchc energy: ', meanNN
        print 'average weight:', meanw
        print 'number of negative energies: ', negative_val

    def Calc_Exchange_Energy(self, crRNA, targetSeq):
        nt_pos = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                  'a': 0, 't': 1, 'c': 2, 'g': 3}
        dG=0
        RNA= ''
        DNA= ''
        for i in range(0, len(crRNA)):
            if i > 0:
                RNA = crRNA[(i-1):(i+1)]
                DNA = targetSeq[(i-1):(i+1)]
                RNA_index = nt_pos[RNA[0]]+4*nt_pos[RNA[1]]
                DNA_index = nt_pos[DNA[0]]+4*nt_pos[DNA[1]]

                dG1 = float(self.decNN[RNA_index][DNA_index])
                if abs(dG1-0.000015) < 1e-6:
                    dG1 = 10000
                    dG1 = 2.3  # during model identification, I set the value of every unknown dG to 0.000015 (if I did not find a value for it)

                pos = 20-i
                w1 = float(self.weights[pos])
                # print 'b1',RNA[0],RNA[1],DNA[0],DNA[1],RNA_index, DNA_index, pos,dG1, w1
            else:
                w1 = 0
                dG1 = 0
            if i < (len(crRNA)-1):
                RNA2 = crRNA[i:(i+2)]
                DNA2 = targetSeq[i:(i+2)]
                RNA_index = nt_pos[RNA2[0]]+4*nt_pos[RNA2[1]]
                DNA_index = nt_pos[DNA2[0]]+4*nt_pos[DNA2[1]]
                dG2 = float(self.decNN[RNA_index][DNA_index])
                if abs(dG2-0.000015) < 1e-6:
                    dG2 = 10000
                    dG2 = 2.3 # during model identification, I set the value of every unknown dG to 0.000015 (if I did not find a value for it)

                pos = 20-i-1
                w2 = float(self.weights[pos])
                # print 'b2',RNA2[0],RNA2[1],DNA2[0],DNA2[1],RNA_index, DNA_index, pos,dG2, w2
            else:
                w2 = 0
                dG2 = 0
            dG += w1*dG1+w2*dG2
        return float(dG)

    def QuickCalc_Exchange_Energy(self,crRNA,TargetSeq):
        nt_pos = {'A': 0, 'T': 1, 'C': 2, 'G': 3,
                  'a': 0, 't': 1, 'c': 2, 'g': 3}
        dG = 0
        RNA = ''
        DNA = ''
        self.nt_mismatch_in_first8 = 0
        for i in range(0, len(crRNA)):
            pos = 20-i
            w1 = self.weights[pos]
            if nt_pos[crRNA[i]] == nt_pos[TargetSeq[i]]:
                dG1 = 0
            else:
                # using a bioinformatics search approach to find sequences with up to x mismatches
                dG1 = 2.3  # kcal/mol
                if pos <= 8:
                    self.nt_mismatch_in_first8 = self.nt_mismatch_in_first8+1
            dG += w1*dG1
        return float(dG)

    def calc_dG_PAM(self, PAM_full_seq):

        # PAM sequence is 5' - N xxx N - 3' where the energy of xxx is listed below.
        # A normal PAM of 'NGG' with 'TC' afterwards would be listed as 'GGT'
        key = PAM_full_seq[1:4]
        if key in self.PAM_energy:
            return self.PAM_energy[key]
        else:
            return 0.0

        # acceptedPAMList=PAM_dic_energy.keys()
        # self.dG_PAM_List=[]
        # self.WarningPAM_List=[]
        # PAMsize=len(self.PAM)
        # for target in self.sequence_list:
            # tPAM=target[-(PAMsize):-1]+target[-1]
            # if tPAM in acceptedPAMList:
                # dGPAM=PAM_dic_energy[tPAM]
                # warning=''
            # else:
                # dGPAM=0
                # warning='N.B'
            # self.dG_PAM_List.append(dGPAM)
            # self.WarningPAM_List.append(warning)

    def calc_dG_exchange(self, guide_sequence, targetSequence):
        self.nt_mismatch_in_first8_list=[]
        if self.quickmode:
            solverfunc=self.QuickCalc_Exchange_Energy
        else:
            solverfunc=self.Calc_Exchange_Energy

        dG_exchange = solverfunc(guide_sequence,targetSequence)

        return dG_exchange

    def calc_dG_supercoiling(self, sigmaInitial, targetSequence):


        sigmaFinal = -0.08
        dG_supercoiling = 10.0 * len(targetSequence) * self.RT * (sigmaFinal**2 - sigmaInitial**2)
        return dG_supercoiling


def get_file_name(args):

    import time
    import datetime
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d_%H_%M')

    return st + "_" + args.target_sequence[:-3] + ".csv"


def get_header(args):

    import time
    import datetime
    ts = time.time()
    timeCurr = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    return [["DATE: ", timeCurr[0:10]],
            ["TIME: ", timeCurr[11:]],
            ["COMMAND:", args],
            ["guide sequence", "position", "target sequence", "dee gee", "partition"]]


def export_file(filedata, filename, args):
    """
    This function write the file data to a csv file
    :param filedata: 
    :param filename: 
    :param args: 
    :return: 
    """
    # TODO: we should use the CSV print lib
    # TODO: we should start the file with a leading title about when this was run

    import csv

    with open(filename, "wb") as csv_output:
        writer = csv.writer(csv_output)
        # writer.writerows(["Run date:",time.time()
        writer.writerows(get_header(args))
        writer.writerows(filedata)

    return


def main():
    import time
    import datetime
    ts = time.time()
    timeCurr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d_%H:%M:%S')

    args = get_options()
    sys.stdout.write("Searching for suitable 'nGG's in: " + '.'.join(args.target_sequence.split('.')[:-1]) + "... ")
    filename = get_file_name(args)
    nggs_list = get_nggs(args.target_sequence)

    Cas9Calculator = clCas9Calculator(args.query, args.model)

    output = []
    for ngg in nggs_list:

        # we need to extract the results of each run into the Output array that we can then print to a file

        sgRNA1 = sgRNA(ngg, Cas9Calculator)
        sgRNA1.run()
        # sgRNA1.printTopTargets()
        # we add this extra line to save the the results to the array
        output.append(sgRNA1.getResults())
        from operator import itemgetter
        output = sortedTargetList = sorted(output, key=itemgetter(4), reverse=True)
        export_file(output, filename, args)

    # sgRNA1.exportAsDill()

    # print sgRNA1

    # PAM='GGA' # NGGA
    # list of all potential on- and off-targets
    # sequence_list=['AGTCCTCATCTCCCTCAAGCCGGA','AGTCCTCATCTCCCTCAAGTCGGA','AGTCCTCATCTCCCTCATGCCGGA']

    # Cas9Calculator=clCas9Calculator(quickmode=True) # using quick approach
    # Cas9Calculator=clCas9Calculator(quickmode=False) # using Invitro or complete model
    # Cas9Calculator=clCas9Calculator(quickmode=False,cModelName='All_dataModel.mat') # using Invitro or complete model
    # Cas9Calculator.loadData(sequence_list,crRNAseq,PAM,True)
    # Cas9Calculator.calcTarget_energy()
    # Cas9Calculator.export_dG() # in an excel file
    # print Cas9Calculator.dG_total_List


main()
