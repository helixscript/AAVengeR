#!/usr/bin/env python3
import os
import argparse
import random
import pandas
import numpy
import re
import sys
import yaml
import subprocess
from Bio.SeqIO import TwoBitIO

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outputDir', type=str, default = '/data/AAVengeR/tests/synData/output', help='Output directory path. Default (/data/AAVengeR/tests/synData/output).', metavar='')
parser.add_argument('-m', '--mode', type=str, default = 'integrase', help='Mode (integrase or AAV). Default (integrase).', metavar='')
parser.add_argument('-n', '--nSites', type=int, default = 100, help='Number of sites to generate. Default (100).', metavar='')
parser.add_argument('-f', '--nFrags', type=int, default = 5, help='Number of fragments for each site (limit 100). Default (5).', metavar='')
parser.add_argument('-p', '--nReadsPerFrag', type=int, default = 25, help='Number of reads per fragment. Default (25).', metavar='')
parser.add_argument('-d', '--distBetweenFragEnds', type=int, default = 25, help='Distance (NT) between fragment break points. Default (25).', metavar='')
parser.add_argument('-s', '--seed', type=int, default = 1, help='Random seed. Default (1).', metavar='')
parser.add_argument('-r', '--refGenomePath', type=str, default = '/data/AAVengeR/data/referenceGenomes/blat/hg38.2bit', help='Path to 2bit reference genome. Default (~/AAVengeR/data/referenceGenomes/blat/hg38.2bit).', metavar='')
parser.add_argument('-g', '--refGenomeID', type=str, default = 'hg38', help='Reference genome ID for sample data table. Default (hg38).', metavar='')
parser.add_argument('-i', '--integraseHMM', type=str, default = 'validation.hmm', help='Name of AAVengeR HMM to use (integrase mode only).  Default (validation.hmm).', metavar='')
parser.add_argument('-t', '--anchorReadStartSeq', type=str, default = 'TCTGCGCGCT', help='Anchor read start sequence filter (AAV mode only). Default (TCTGCGCGCT).', metavar='')
parser.add_argument('-x', '--R1_length', type=int, default = 150, help='Total length of R1 reads. Default (150).', metavar='')
parser.add_argument('-y', '--R2_length', type=int, default = 150, help='Total length of R2 reads. Default (150).', metavar='')
parser.add_argument('-z', '--I1_length', type=int, default = 12, help='Total length of I1 reads. Default (12).', metavar='')
parser.add_argument('-e', '--percentGenomicError', type=float, default = 0, help='Percent gDNA error (0.0 - 1.0) to simulate in R1 and R2 reads. Default (0).', metavar='')
parser.add_argument('-c', '--positionChatterSD', type=float, default = 0.50, help='StdDev of Gaussian centered on expected fragment ends used to simulate position chatter. (Default 0.50).', metavar='')
parser.add_argument("-u", "--softwareDir", type=str,    default = '/data/AAVengeR',        help = "Full path to AAVengeR software as seen from Docker.", metavar = "")
parser.add_argument("-a", "--threads", type =int,   default = 20,            help = "Number of computation threads to use.", metavar = "")
parser.add_argument("-k", "--alignmentMinPercentID", type=float,  default = 95,     help = "BLAT percent seq id.", metavar = "")
parser.add_argument("-w", "--alignmentBlatRepMatch", type=float,  default = 5000,   help = "BLAT repMatch parameter value.", metavar = "")

args = parser.parse_args()
args.refGenomePath = os.path.expanduser(args.refGenomePath)

# Test inputs 
if not os.path.exists(args.refGenomePath):
  print('Error - refGenomePath does not')
  sys.exit(1)

if args.mode not in ['integrase', 'AAV']:
  print('Error - mode must be set to "integrase" or "AAV".')
  sys.exit(1)

if args.nFrags > 50:
  print('Error - nFrags must be set to a value no more than 50.')
  sys.exit(1)

if args.percentGenomicError < 0 or args.percentGenomicError > 1:
  print('Error - percentGenomicError must be set to a value between 0 and 1.')
  sys.exit(1)

if args.positionChatterSD < 0 or args.positionChatterSD > 5:
  print('Error - positionChatterSD must be set to a value between 0 and 5')
  sys.exit(1)

if(args.nFrags * args.distBetweenFragEnds > 500):
  print('Error - the number of requested fragments x the distance between fragment ends can not exceed 500.')
  sys.exit(1)

# Set random seed.
random.seed(args.seed)

# Create output directory.
if not os.path.isdir(args.outputDir):
  os.mkdir(args.outputDir)

if not os.path.exists(args.outputDir):
  print('Error - could not create output directory.')
  sys.exit(1)

# Empty out output dir if it contains a previous result.
files = [os.path.join(args.outputDir, 'R1.fastq.gz'),
         os.path.join(args.outputDir, 'R2.fastq.gz'),
         os.path.join(args.outputDir, 'I1.fastq.gz'),
         os.path.join(args.outputDir, 'testVector.fasta'),
         os.path.join(args.outputDir, 'sampleData.tsv'),
         os.path.join(args.outputDir, 'config.yml'),
         os.path.join(args.outputDir, 'truth.tsv')]

for f in files:
  if os.path.exists(f):
    os.remove(f)


# Linkers to be used for faux sites.
remnant0 = 'TCTGCGCGCTCGCTCGCTCA' # To be used for integrase mode.
remnant  = 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTG' 
linker   = 'GAACGAGCACTAGTAAGCCCNNNNNNNNNNNNCTCCGCTTAAGGGACT' 


# Helper functions.
def flatten(xss):
    return [x for xs in xss for x in xs]

def randomBarCode(length = 12):
  NTs = ['A', 'T', 'C', 'G']
  return(''.join(random.choices(NTs, k = length)))

def compSeq(seq):
  compBases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
  return(''.join([compBases[x] for x in seq]))

def revCompSeq(seq):
  return(compSeq(seq)[::-1])

def buildGaussDistSample(pos, n, popSize = 1000, sd = 1):
  nums = []
  for x in range(popSize):
    nums.append(round(random.gauss(pos, sd)))
  return(random.sample(nums, n))

def mutateAtPos(string, index):
    NTs = ['A', 'T', 'C', 'G']
    exclude = list(string[index])
    NTs =[x for x in NTs if x not in exclude]
    return string[:index] + random.sample(NTs, 1)[0] + string[index + 1:]

def simulateError(seq):
  nPos = round(args.percentGenomicError * len(seq))
  orgSeq = seq

  locs = random.sample(list(range(1, len(seq))), nPos)
  for x in locs:
    seq = mutateAtPos(seq, x)

#  ed = Levenshtein.distance(orgSeq, seq)
#  if ed >= 5:
#      print('nPos: ', nPos, ' orgLen: ', len(seq), ' new length: ', len(seq), ' Edit dist: ', ed)

  return(seq)


# Build a collection of scrambled remnants to draw from. 
# Pieces will be between 12 and 36 NTs with a 50% of being flipped to RC.
start = []
end = []

print('Building remnant library.')

for i in range(500):
  while True:
    a = random.randint(12, len(remnant)-12)
    b = random.randint(a+12, len(remnant))
    
    if (b - a + 1) >= 12 & (b - a + 1) <= 36:
      start.append(a)
      end.append(b)
      break

pieces = []
for i in range(len(start)):
  pieces.append(remnant[start[i-1]:end[i-1]])

remnants = []
for i in range(500):
   n = random.randint(0, 3)

   if n == 0:
     remnants.append(remnant[0:random.randint(12, 24)])
   else:
      # Start remnant with a variable length fragment starting from the beginning.
      o = [remnant[0:random.randint(12, 24)]]

      # Randomly draw 1 - 3 remnant pieces.
      k = random.sample(range(1, len(pieces)), n)
      p = [pieces[x] for x in k]
    
      # Randomly flip pieces to RC.
      p2 = []
      for r in p:
         if random.randint(0, 1) == 0:
            p2.append(r)
         else:
            p2.append(revCompSeq(r))

      o.append(p2)
      o = flatten(o)
      remnants.append(''.join(o))

if args.mode == 'integrase':
   remnants = [remnant0]


# Read in genomic data pointers.
print('Reading reference genome data.')
_handle = open(args.refGenomePath, 'rb')          # keep file open for lazy-loading
g = TwoBitIO.TwoBitIterator(_handle)

# Create a list of allowed chromosomes from which to sample.
o = list(range(1, 100))
o = [str(x) for x in o]
o.extend(['X'])  # Exclude Y, often full of repeats and poorly characterized.
allowedChromosomes = ['chr' + x for x in o]
chromosomes = []
for c in g.keys():
    if c in allowedChromosomes:
        chromosomes.append(c)

# Read select chromosomes into memory, slow.
d = {}
for c in g.keys():
    if c in chromosomes:
        print('  Reading chromosome', c, 'into memory.')
        d[c] = str(g[c].seq)   # g[c] is a SeqRecord; use .seq

del g
_handle.close()

# Select chromosomes for sites.
chromosomes = random.choices(chromosomes, k = args.nSites)

# Compile data to build sites.
class site:
    def __init__(self, chr, pos, seq, strand):
        self.chr = chr
        self.pos = pos
        self.seq = seq
        self.strand = strand

sites = []
print("Compiling site building data.")
for c in chromosomes:   
  # Select a position within the chromosome excluding the ends.
  pos = random.randint(1+10000, len(d[c])-10000)
  seq = 'N'
    
  # Select a 1000 NT chunk from the genome that does not contain any Ns.  
  # This R1 and R2 reads will be created from the both sides of this large chunk.
  while True:
    seq = d[c][pos:pos+1000]
    if 'N' not in seq:
      break
      
    # Select new position if previous yielded a sequence with Ns.
    pos = random.randint(1+10000, len(d[c])-10000)

  strand = random.sample(['+', '-'], 1)[0]

  if strand == '-':
    pos = pos + 1000
    seq = revCompSeq(seq)

  sites.append(site(c, pos, seq, strand))


# Loop through site objects and build read sequences.
print("Building site reads.")
for s in sites:
  s.readIDs = []
  s.UMIs = []
  s.fragmentSeqs = []
  s.R1 = []
  s.R2 = []

  nReads = args.nFrags * args.nReadsPerFrag

  # Start within excised gDNA chunks to allow chatter around start position.
  # Create start positions around position 10.
  startPositions = buildGaussDistSample(10, nReads, sd = args.positionChatterSD)

  # Create end positions simulating sonic breaks.
  baseFragWidth = 500
  readNum = 1

  # Create alt start postion for read ids to correct for zero-based 2bit system and starting 10 NT within seq blocks.
  # Integrase vs AAV mode invoke different gDNA duplication corrections in buildSites.
  if s.strand == '+':
    if(args.mode == 'integrase'):
       s.pos2 = s.pos + 12
    else:
       s.pos2 = s.pos + 10
  else:
    if(args.mode == 'integrase'):
      s.pos2 = s.pos - 11
    else:
      s.pos2 = s.pos - 9

  endPositions = []

  for x in range(1, args.nFrags + 1):
    # Define a frag specific UMI sequence.
    UMI = ''.join(random.choices(['A', 'T', 'C', 'G'], k = 12))

    for i in range(1, args.nReadsPerFrag + 1):
      s.readIDs.append(s.chr + s.strand + str(s.pos2) + '_read' + str(readNum) + '_frag' + str(x))
      s.UMIs.append(UMI)
      readNum += 1

    offSet = buildGaussDistSample(baseFragWidth + (x * args.distBetweenFragEnds), args.nReadsPerFrag, sd = args.positionChatterSD)
    endPositions.append(offSet)

  endPositions = flatten(endPositions)

  # Extract substrings for fragment sequences.
  for x in range(len(s.readIDs)):
     s.fragmentSeqs.append(s.seq[startPositions[x]-1:endPositions[x]-1])

  # Randomly select a remnant. The same remnant will always be selected when in 'integrase' mode.
  s.remnant = remnants[random.randint(0, len(remnants)-1)]

  # Build R2 and R1 read sequences.
  for x in range(len(s.fragmentSeqs)):
    gDNA = s.fragmentSeqs[x][0:args.R2_length - len(s.remnant)]

    if(args.percentGenomicError > 0):
      saved_state = random.getstate()
      gDNA = simulateError(gDNA)
      random.setstate(saved_state)

    s.R2.append(s.remnant + gDNA)

  for x in range(len(s.fragmentSeqs)):
    gDNA = revCompSeq(s.fragmentSeqs[x]) [0:args.R1_length - len(linker)]

    if(args.percentGenomicError > 0):
      saved_state = random.getstate()
      gDNA = simulateError(gDNA)
      random.setstate(saved_state)

    s.R1.append(linker + gDNA)


# Build sample data file.
print("Building sample table.")
sampleData = pandas.DataFrame({
  'replicate': [1,2,3,4] * 9,
  'subject': flatten([['subjectA'] * 12, ['subjectB'] * 12, ['subjectC'] * 12]),
  'sample': flatten([['sample1'] * 4, ['sample2'] * 4, ['sample3'] * 4]*3),
  'trial': 'test',
  'index1Seq': [randomBarCode() for i in range(36)],
  'adriftReadLinkerSeq': linker,
  'refGenome': args.refGenomeID
})


# Add mode specific columns.
if(args.mode == 'integrase'):
  sampleData['leaderSeqHMM'] = args.integraseHMM
  sampleData['flags'] = 'IN_u5'
  sampleData['vectorFastaFile'] = 'none.fasta'
else:
  sampleData['anchorReadStartSeq'] = args.anchorReadStartSeq
  sampleData['flags'] = 'AAV'
  sampleData['vectorFastaFile'] = 'validationVector.fasta'

# Split sample table and sites into sample groups.
sampleGroups = sampleData.groupby(['trial', 'subject', 'sample'])
siteGroups = numpy.array_split(sites, sampleGroups.ngroups)

z = 0
truths = []

# Create a dictionary of fragment ids to replicate numbrtd (1 - 4). 
# (!) Only goes out to frag100 - need to limit user frag requests. 
fragToRep = dict(zip(['frag' + str(n) for n in list(range(1,101))], list(range(1, 5))*20))

print("Building I1 reads and writing all reads to output.")
for (key, group) in sampleGroups:
    # Create list of replicate bar codes associated with this sample.
    barCodes = list(group['index1Seq'])
    
    # Loop through sites in c
    for x in siteGroups[z]:

      # Keep fragments within replicates since they are standardized within replicates.
      fragIDs = [re.search('frag\\d+', m).group(0) for m in x.readIDs]
      reps = [fragToRep[n] for n in fragIDs]
      barCodesToUse = [barCodes[x-1] for x in reps]

      # Add entry to truths.
      truths.append({'trial': key[0], 'subject': key[1], 'sample': key[2], 'posid': x.readIDs[i].split('_')[0], 
                    'nReads': args.nFrags * args.nReadsPerFrag, 'nFrags': args.nFrags, 'nUMIs': args.nFrags, 
                    'leaderSeq': x.remnant})
      
      for i in range(len(x.readIDs)):
        # Replace NNN in linker with fragment UMIs, build R1 FASTQ, write.
        x.R1[i] = x.R1[i].replace('NNNNNNNNNNNN', x.UMIs[i])
        o = ['@' + x.readIDs[i], x.R1[i], '+', ''.join('?'*len(x.R1[i]))]

        with open(os.path.join(args.outputDir, 'R1.fastq'), 'a') as f:
           for line in o:
             f.write("%s\n" % line)

        # Build and write R2.
        o = ['@' + x.readIDs[i], x.R2[i], '+', ''.join('?'*len(x.R2[i]))]
        
        with open(os.path.join(args.outputDir, 'R2.fastq'), 'a') as f:
         for line in o:
           f.write("%s\n" % line)

        # Build and write I1.
        # Select a barcode from the ones associated with this sample.
        b = barCodesToUse[i]
        o = ['@' + x.readIDs[i], b, '+', ''.join('?'*len(b))]

        with open(os.path.join(args.outputDir, 'I1.fastq'), 'a') as f:
         for line in o:
           f.write("%s\n" % line)
    z += 1     

print("Writing tabular outputs.")
# Write out sample table.
sampleData.to_csv(os.path.join(args.outputDir, 'sampleData.tsv'), sep='\t', index=False) 

# Write out truth table.
t = pandas.DataFrame.from_records([truth for truth in truths])
t.to_csv(os.path.join(args.outputDir, 'truth.tsv'), sep='\t', index=False) 

# Compress fastq files.
os.system('gzip ' + args.outputDir + '/*.fastq')

print('Preparing config file.')

config_path = os.path.join(args.softwareDir, "config.yml")
with open(config_path, "r", encoding="utf-8") as f:
    configData = yaml.safe_load(f)

configData['mode'] = args.mode
configData['outputDir'] =  os.path.join(args.outputDir, 'output')
configData['softwareDir'] = args.softwareDir
configData['core_CPUs'] = args.threads
configData['database_samplesConfigGroup'] = 'none'
configData['demultiplex_anchorReadsFile'] = os.path.join(args.outputDir, 'R2.fastq.gz')
configData['demultiplex_adriftReadsFile'] = os.path.join(args.outputDir, 'R1.fastq.gz')
configData['demultiplex_index1ReadsFile'] = os.path.join(args.outputDir, 'I1.fastq.gz')
configData['demultiplex_sampleDataFile'] = os.path.join(args.outputDir, 'sampleData.tsv')
configData['alignReads_genomeAlignment_minPercentID'] = args.alignmentMinPercentID
configData['alignReads_genomeAlignment_blatRepMatch'] = args.alignmentBlatRepMatch
configData['modules'] = ['core']

# Set these valuesto 1 to prevent errors parseing values such as 1,2,3. These paramters are not used for core testing.
configData['anchorReadStartSeqs_windows'] = 1
configData['anchorReadRearrangements_seeds'] = 1
configData['anchorReadRearrangements_rarifactionLevels'] = 1
configData['pullDatabaseFragments_trialSubjectSamples'] = 1

with open(os.path.join(args.outputDir, "config.yml"), "w", encoding="utf-8") as f:
    yaml.safe_dump(
        configData, f,
        sort_keys = False,           # keep your key order
        default_flow_style = False,  # block style
        allow_unicode = True
    )

print('Running AAVengeR ...')

comm = ['docker', 'run', '--rm', '--mount', 'type=bind,source=/data,target=/data', '-e', 'CONFIG_PATH=' + os.path.join(args.outputDir, 'config.yml'), 'aavenger']

result = subprocess.run(comm, check=True)

print('Running evalSynDataResult ...')

comm = [os.path.join(args.softwareDir, 'tests', 'synData', 'evalSynDataResult.R'),
        '--outputDir',    os.path.join(args.outputDir, 'eval'),
        '--sitesFile',    os.path.join(args.outputDir, 'output', 'core', 'sites.rds'),
        '--multiHitFile', os.path.join(args.outputDir,  'output', 'core', 'multiHitClusters.rds'),
        '--truthFile',    os.path.join(args.outputDir, 'truth.tsv'),
        '--siteWidth', '3' ]

result = subprocess.run(comm, check=True)

print('Done.')
