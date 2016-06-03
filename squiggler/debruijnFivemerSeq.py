import sys, itertools
sys.path.append("/Users/johnurban/MOOCs/Rosalind/Stepic/bioinfoPevznerCompaeu/molecularClock/")
sys.path.append("/Users/johnurban/MOOCs/Rosalind/Stepic/bioinfoPevznerCompaeu/genomeAssembly/")
sys.path.append("/Users/johnurban/MOOCs/Rosalind/Stepic/bioinfoPevznerCompaeu/originFinding/")
from molClock import *
from genomeAssemblyTools import *
from oriSeekingTools import *

fivemers = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT")]
g=deBruijnBuildFromKmers(fivemers)
e=g.getEdgeList()
c = findEulerianCycle(e)
s = getSequenceFromEulPath(c)

##failct=0
##for fivemer in fivemers:
##    if fivemer not in s:
##        print fivermer, "not in sequence"
##        failct+=1
##        #this should never print a fivemer, it is just a check
##print failct

        
print s

print
print "For 250 bp oligos overlapping by 50 bp, get..."
s250 = [s[:250], s[200:450], s[400:650], s[600:850], s[800:]]
print ("\n\n").join(s250)
print
print 
print "...and their reverse complements to anneal in tube."
print
print ("\n\n").join([reverseComplement(seq) for seq in s250])

print
print
print "Alt. get staggered overlapping oligos that can be ligated together to form the 1028 bp sequence"
print "---------------- ---------------- -----------------"
print "--------- ------------- ---------------- ----------\n"
print "ligate"
print "----------------*----------------*-----------------"
print "---------*-------------*----------------*----------\n"
print "Walaaa"
print "---------------------------------------------------"
print "---------------------------------------------------"
    
seq='TTTTATTTAATTTCTTTTCATTTGTTTTGATTATACTATAGTATATTAAACTAAAGTAAATTACCCTACCGTACCTTACGCTACGGTACGTTACTCTACTGTACTTTACAATACACTACAGTACATTAGCCTAGCGTAGCTTAGGCTAGGGTAGGTTAGTCTAGTGTAGTTTAGAATAGACTAGAGTAGATCCTATCCAACCCAAGCCAATCCCCCGCCCCTCCCGGCCCGTCCCTGCCCTTCCCACCCCAGCCCATACCATCCGCGCCGCTCCGGGCCGGTCCGTGCCGTTCCGAACCGACCCGAGCCGATACGATCGTATCGAAGCGAATCGCCTCGCGGCGCGTCGCTGCGCTTCGCAACGCACCGCAGCGCATAGCATCGGCTCGGGGCGGGTCGGTGCGGTTCGGAACGGACCGGAGCGGATAGGATCTTATCTATGCTATTCTAACCTAAGCTAATATAATCTCCTCTCGTCTCTGCTCTTCTCAACTCACCTCAGCTCATATCATCTGGCCTGGGCTGGTATGGTCTGTCCTGTGCTGTTATGTTCTGAACTGACCTGAGCTGATATGATCAAGGCAAGTCAATGCAATTCAAAAACAAACCAAAGAAAAGCAAATAAAATCACGCCACGGCACGTAACGTCACTCCACTGCACTTAACTTCACAACACACCACAGAACAGCACATAACATCAGGCCAGGGAAGGGCAGGTAAGGTCAGTCCAGTGAAGTGCAGTTAAGTTCAGACAAGACCAGAGAAGAGCAGATAAGATGCCTTGCCATGCGTGGCGTTGCGACGCGAGGCGATGGAATGGCTTGGCATGGGGGTGGGTTGGGACGGGAGGGGATGTATTGTAATGTCGTGTCTTGTCATGTGGTGTGTTGTGACGTGAGGTGATGAATTGAAACGAAAGGAAATGACTCGACTGGACTTGACACGACAGGACATGAGTCGAGTGGAGTTGAGACGAGAGGAGATTCCTTTCCATTCGTTTCGATTGCTTTGCATTGGTTTGGATTTTT'
