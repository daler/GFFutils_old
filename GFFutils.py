"""
Module for complex interaction with a GFF file.

See http://github.com/daler/GFFutils for source and documentation.
"""
import os
import sqlite3
import itertools
import sys
import time
import tempfile
import mmap
import string
import copy
import gzip

class GFFFeature(object):
    """
    Class to represent a GFF feature (gene, mRNA, exon, etc)
    """
    class Attributes(object):
        '''
        Simple wrapper class to provide easy access to data from GFF
        "attribute" fields.
        '''
        def __init__(self):
            self._attrs = []  # will hold a list of attributes added to the object.

    def __init__(self, chr, source, featuretype, start, stop,
                 value,strand,phase,attributes,strvals=False):
        """
        *chr*
        
            Chromosome
        
        *source*
        
            Source of the data
        
        *featuretype*
        
            type of the feature ('gene', 'mRNA', 'CDS', etc)
        
        *start*
        
            Start position on the chromosome
        
        *stop*
            
            Stop position on the chromosome
        
        *value*
        
            Value for the feature
        
        *strand*
        
            Strand of the feature, '+' or '-'
        
        *phase*
        
            Phase of the feature if it's a CDS
        
        *attributes*
        
            A semicolon-delimited list of "field=data" strings.  For example, 
            'ID=FBtr0000123;Parent=FGgn0001;Name=transcript 1'
        
        *strvals*
        
            By default, GFFFeature objects will have their attributes typecast
            to integers, floats, etc.  However if *strvals* is True, ALL
            attributes will be strings.  For example, with *strvals=False*,
            GFFFeature.start and GFFFeature.end be integers (useful for
            downstream work like ``featurelen = feature.stop - feature.start``)
            but if *strvals=True* they will be strings. 
        
            Setting *strvals=True* will speed up parsing.
        """

        if not strvals: # do typecasting
            self.chr=chr
            self.source=source
            self.featuretype=featuretype
            try:
                self.start=int(start)
            except ValueError:
                self.start = None
            try:
                self.stop=int(stop)
            except ValueError:
                self.stop = None
            if value is not None:
                try:
                    self.value=float(value)
                except ValueError,TypeError:
                    self.value=None
            else:
                self.value = None
            self.strand=strand
            if phase is not None:
                try: 
                    self.phase=int(phase)
                except ValueError,TypeError:
                    self.phase=None
            else:
                self.phase = None
          
        if strvals: # no typecasting, save everything as a string.
            self.chr=chr
            self.source=source
            self.featuretype=featuretype
            self.start=start
            self.stop=stop
            self.value=value
            self.strand=strand
            self.phase=phase

        self._parse_attributes(attributes)

    def _parse_attributes(self,attributes):  
        self._strattributes = attributes # keep track of these for later printing out.
        if attributes is None:
            self._strattributes = ''
        # parse "attributes" field of the GFF line and insert them into an Attributes 
        # object.
        self.attributes = GFFFeature.Attributes()
        if attributes is not None:
            items = attributes.split(';')
            for item in items:
                if len(item) > 0:
                    field,value = item.split('=')
                field = field.strip()
                values = value.split(',')
                values = [i.strip() for i in values]
                setattr(self.attributes,field,values)

                # Keep track inside the Attributes object of what you added to it
                self.attributes._attrs.append(field)
        


    @property
    def id(self):
        try:
            return self.attributes.ID[0]
        except AttributeError:
            return None

    def add_attribute(self,attribute,value):
        """
        Add an attribute to this feature.

        *value* can be of any type; if it can't be sliced (and it's not a
        string) then it's converted into a list automatically.  
        """
        try:
            value[0]
        except TypeError:
            # unsubscriptable
            value = [value]

        # Strings are subscriptable, but should be wrapped as a list.
        if type(value) == str:
            value = [value]
        setattr(self.attributes,attribute,value)
        self.attributes._attrs.append(attribute)
        attr = ','.join(value)
        self._strattributes+=';%s=%s' % (attribute,attr)

    def remove_attribute(self,attribute):
        """
        Delete attribute from this feature.
        """
        delattr(self.attributes, attribute)
        self.attributes._attrs.remove(attribute)
        ind1 = self._strattributes.find(attribute)
        ind2 = self._strattributes.find(';',ind1)
        self._strattributes = self._strattributes[:ind1]+self._strattributes[ind2:-1]

    def to_bed(self,fieldcount=3):
        """
        Returns the feature as a BED format line, with number of fields
        *fieldcount*.  Default is 3, which is chrom, start, stop.  Up to BED-6
        is supported.
        """
        attrs = ['chr','start','stop','id','value','strand']
        fields = []
        for i in range(fieldcount):
            fields.append(str(getattr(self,attrs[i])))
        return '\t'.join(fields)+'\n'

    def __repr__(self):
        return "%s %s '%s': %s:%s-%s (%s)" % (self.__class__.__name__,self.featuretype,self.id,self.chr,self.start,self.stop, self.strand)

    def __len__(self):
        return self.stop-self.start
   
    @property
    def TSS(self):
        """The transcription start site of the feature.  This is simply the
        strand-specific start of the feature.  Returns None if no strand
        specified.
        """
        if self.strand == '+':
            return self.start
        if self.strand == '-':
            return self.stop
    
    @property
    def midpoint(self):
        """Convenience accessor for getting the midpoint of a feature.
        """
        return self.start + (self.stop-self.start)/2

    def tostring(self):
        """Prints the GFF record suitable for writing to file (newline included).
        
        Since the string output is reconstructed based on the current contents
        of the GFFFeature, (attributes are reconstructed as well), this
        provides an easy means of editing GFF files, e.g.::
            
            # shorten all CDS features by 10 bp, rename features to "*.short",
            # and write only these shortened CDSs to file.

            fout = open('out.gff')
            for feature in GFFFile('in.gff'):
                if feature.featuretype != 'CDS':
                    continue
                feature.stop -= 10
                feature.attributes.ID += '.short'
                fout.write(feature.tostring())
            fout.close()

        In the interest of speed, does not do error-checking.
        """
        # Reconstruct the attributes field
        attributes = []
        for attr in self.attributes._attrs:
            values = getattr(self.attributes,attr)
            values = map(str,values)
            values = ','.join(values)
            attributes.append(attr+'='+values)
        attributes = ';'.join(attributes)

        items = [self.chr, 
                 self.source,
                 self.featuretype,
                 self.start, 
                 self.stop, 
                 self.value, 
                 self.strand, 
                 self.phase,
                 attributes]

        printables = []
        for item in items:
            if item is None:
                printables.append('.')
            else:
                printables.append(str(item))
        return '\t'.join(printables).rstrip()+'\n'

class GFFFile(object):
    """Iterator object, that moves through features in a GFF-format file.  A
    new gfffeature object is created for each line.
    
    Usage::

        for feature in gfffile('a.bed'):
            print feature.chr
            print feature.start
            print feature.stop
            print feature.featuretype
            print feature.desc
            print 'length:', feature.stop-feature.start
        """

    featureclass = GFFFeature

    def __init__(self,f,strvals=False):
        """
        *f*

            GFF filename.  Can be a string filename (.gz files detected via extension)
            or, if *f* is not a string type, it can be a file-like object (sys.stdout, 
            already-open file, StringIO, etc).

        *strvals*

            By default, GFFFeature objects will have their attributes typecast
            to integers, floats, etc.  If *strvals* is True, ALL attributes
            will be strings.  For example, with *strvals=False*,
            GFFFeature.start and GFFFeature.end be integers (useful for
            downstream work) but if *strvals=True* they will be strings. 

            Setting *strvals=True* will speed up parsing.
        """
        if type(f) is str:
            self.stringfn = True
            if os.path.splitext(f)[-1] == '.gz':
                self.file = gzip.open(f)
            else:
                self.file = open(f)
        else:
            self.stringfn = False
            self.file = f
        self.strvals = strvals

    def __iter__(self):
        """Iterator function.  Yields a gfffeature object for each line."""
        f = self.file
        for line in f:
            # You've reached the end of the GFF file; this represents the start of 
            # the optional sequence section
            if line.startswith('>'):
                raise StopIteration
            line = line.rstrip()
            if line.startswith('#') or len(line) == 0:
                continue
            L = line.rstrip().split('\t')
            args = [None for i in range(9)]
            args[:len(L)] = L
            args.append(self.strvals)
            if self.__class__.featureclass == GFFFeature:
                yield self.__class__.featureclass(*args)
            if self.__class__.featureclass == GTFFeature:
                args.insert(0,None)
                yield self.__class__.featureclass(*args)

        # close up shop when done.
        if self.stringfn:
            f.close()
    
    def __eq__(self,other):
        if self._strattributes != other._strattributes:
            return False
        for attr in ['id','featuretype','start','stop','chr']:
            if getattr(self,attr) != getattr(other,attr):
                return False
        return True

    def __repr__(self):
        return 'gfffile object (file=%s)' % (self.file)


exon_numbers = {}

class GTFFeature(GFFFeature):
    """
    Class to represent a GTF feature and its annotations. Subclassed from GFFFeature.
    """
    def __init__(self, id, chr, source, featuretype, start, stop,
                 value,strand,phase,attributes,strvals=False):
        #print id,chr,source,featuretype,start,stop,value,strand,phase,attributes
        GFFFeature.__init__(self,chr,source,featuretype,start,stop,value,strand,phase,attributes,strvals)
        self.add_attribute('ID',id)
    
    def tostring(self):
        """
        Prints the GTF record suitable for writing to file (newline
        included).
        
        In the interest of speed, does not do error-checking.

        AB000123    Twinscan     CDS    193817    194022    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2"; 
        """
        # Genes and transcripts are not explicitly written in GTF files.
        if self.featuretype == 'gene':
            return ''
        if self.featuretype == 'mRNA':
            return ''
        
        # Reconstruct the attributes field
        attributes = ''
        for attr in self.attributes._attrs:
            values = getattr(self.attributes,attr)
            if type(values) is list:
                values = ','.join(map(str,values))
            attributes += attr+' '+'"'+str(values)+'"; '

        items = [self.chr, 
                 self.source,
                 self.featuretype,
                 self.start, 
                 self.stop, 
                 self.value, 
                 self.strand, 
                 self.phase,
                 attributes]

        printables = []
        for item in items:
            if item is None:
                printables.append('.')
            else:
                printables.append(str(item))
        return '\t'.join(printables).rstrip()+'\n'

    def add_attribute(self,attribute,value):
        """
        Add an attribute to this feature.

        *value* can be of any type; if it can't be sliced (and it's not a
        string) then it's converted into a list automatically.  
        """
        try:
            value[0]
        except TypeError:
            # unsubscriptable
            value = [value]

        # Strings are subscriptable, but should be wrapped as a list.
        if type(value) == str:
            value = [value]
        setattr(self.attributes,attribute,value)
        self.attributes._attrs.append(attribute)
        assert len(value) == 1, 'Multiple values for one attribute not currently supported for GTF features'
        self._strattributes+=';%s "%s"' % (attribute,value[0])
    

    def remove_attribute(self,attribute):
        """
        Delete attribute from this feature.
        """
        delattr(self.attributes, attribute)
        self.attributes._attrs.remove(attribute)
        ind1 = self._strattributes.find(attribute)
        ind2 = self._strattributes.find(';',ind1)
        self._strattributes = self._strattributes[:ind1]+self._strattributes[ind2:]

    def _parse_attributes(self,attributes):
        """
        Parse the attributes.  This is where GTF differs from GFF format.
        """
        self._strattributes = attributes # keep track of these for later printing out.

        # parse "attributes" field of the GTF line and insert them into an Attributes 
        # object.
        self.attributes = GTFFeature.Attributes()
        if attributes is not None:
            items = attributes.split(';')
            for item in items:
                if len(item) == 0:
                    continue
                field,value = item.strip().split()
                value = value.replace('"','')
                try:
                    value = float(value)
                except:
                    pass
                setattr(self.attributes,field,value)

                # Keep track inside the Attributes object of what you added to it
                self.attributes._attrs.append(field)
        
           

class GTFFile(GFFFile):
    """Iterator object, that moves through features in a GTF-format file.  A
    new GTFFeature object is created for each line.  Subclassed from GFFFile.
    
    Usage::

        for feature in GTFFile('a.bed'):
            print feature.chr
            print feature.start
            print feature.stop
            print feature.featuretype
            print feature.desc
            print 'length:', feature.stop-feature.start
        """
    featureclass = GTFFeature

    def __repr__(self):
        return 'GTFFile object (file=%s)' % (self.file)

def gtf2gff3(gtf, gff):
    """
    Calls the perl script.
    """
    perlscript = os.path.join(exe_path,'..','scripts','gtf2gff3.pl')
    perlscript = os.path.abspath(perlscript)
    cmd = [perlscript, gtf,'>', gff]
    p = subprocess.call(cmd).communicate()
    #os.system(' '.join(cmd))

class Genome:
    """
    Wrapper class for quickly getting a sequence within a chromosome. Inspired by ERANGE.

    Creates a memory-map of the fasta file.  It indexes where the newlines are
    in the file and records how many bits it is into the file.  Access is then
    pretty quick, since we can seek to that bit in the file rather quickly.
    
    Example usage::

        >>> g = Genome('dm3.fa')
        >>> nucleotides = g.sequence('chr2L',12000,13000)
    
    This implementation does NOT depend on BioPython, but converting to a
    BioPython Seq object is pretty easy:

        >>> from Bio.Seq import Seq
        >>> seq = Seq(nucleotides)
        >>> revcomp = seq.reverse_complement()
        >>> protein = seq.translate()

    """
    def __init__(self,fn):
        """
        *fn*

            fasta file with sequences each on a single line.
        """
        # Assume file has a single newline, and it's what separates the
        # description from the sequence.  
        #
        # e.g., 
        #
        # >chr2L
        # AAAGATCTGACTGACCGCGCGGGATATCGCGCATGCTAC...
        # > chr2R
        # ACGATCGCGCGCAATAATTTATATGCGACTAGCTGTAGC...

        # Keep track of the starting position of each chromosome in self.startinds
        self.startinds = {}
        self.chromfiles = {}
        f = open(fn)
        m = mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)
        ind1 = 0
        m.seek(0)
        while True:
            ind1 = m.find('>')
            if ind1 == -1:
                break
            m.seek(ind1)
            ind2 = m.find('\n')
            print ind1
            print ind2
            chrom = m[ind1+1:ind2]
            m.seek(ind2+1)
            print chrom
            self.startinds[chrom] = ind2+1
        self.mmap = m

    def sequence(self,chrom,start,stop):
        """
        Returns the sequence for the position requested.
        """
        i = self.startinds[chrom]
        start = i + start
        length = i + stop - start
        self.mmap.seek(start)
        seq = self.mmap.read(length)

        # Since we're reading right from the file, seq may have newlines in it.
        # So count the number of newlines, get another, longer sequence (longer
        # by the number of newlines) and return the new longer sequence with
        # the newlines stripped.
        newlines = seq.count('\n')
        if newlines == 0:
            return seq
        
        if seq.startswith('\n'):
            # then shift it back one
            start -= 1
        
        self.mmap.seek(start)
        seq = self.mmap.read(length+newlines)
        return seq.replace('\n','')

#TODO . . . stopped working here.
def splicejunctions(genomefasta, gffdbfn, splicewidth, gene=None, transcript=None):
    """
    Returns an iterator of dictionaries containing splice junction sequences. 
    
    Width of splice junctions is determined by *splicewidth*. 

    Needs as input a FASTA file of the genome (*genomefasta*) which will be memory-mapped, 
    and a GFF database that GFFDB can use (*gffdbfn*)

    If an ID for *gene* is provided, only that gene's transcripts' splice junctions will be returned.
    If an ID for *transcript* is provided, only that transcript's splice junctions will be returned.
    
    Each dictionary returned has two keys, 'gene' and 'transcript'.  'gene' contains the gene name, and
    'transcripts' contains a dictionary of transcripts, each of which contain a list of splice junction 
    sequences::

        d = {gene='FBgnNNNN',
             transcripts = {FBtrNNNN1=['ATCGCGTTGC','ACGCGCGCGT', 'CCCCTTTTAACT'],
                            FBtrNNNN2=['CTCGCGCTCC','CCCCACACAA', 'CGCCGTCTGACC'],
            }

    Splices are reverse-complemented if on the minus strand.
    """

    G = GFFDB(gffdbfn)
    genome = Genome(genomefasta)
    
    d = {}
    if gene is not None:
        child_transcripts = G.children(gene)

        for child in child_transcripts:
            if child.featuretype != 'mRNA':
                continue
            exons = [i for i in G.children(child) if i.featuretype=='exon']
            
            #exons should already be sorted upon return from database

            for i in exons:
                for j in exons:
                    part1_start = i.stop - halfwidth
                    part1_stop = i.stop
                    part1_seq = genome.sequence(chrom,part1_start,part1_stop)
                    part2_start = j.start
                    part2_stop = j.start + halfwidth
                    part2_seq = genome.sequence(chrom,part2_start,part2_stop)
                    splice_seq = part1_seq + part2_seq
                    
            

def splices(seqfile,gffdbfn,readlen,fout):
    """
    Currently writes splices out to a file.  Would be cool if you could iterate
    through splices from DB.

    seqfile is a fasta-formatted file with the original chromosome sequences.

    gffdbfn is the filename of a gff database populated by gffdb()

    readlen is the length of reads

    fout can either be a string or a file-like object.
    """
    separator = ':'
    close_fout = False
    if type(fout) is str:
        close_fout = True
        fout = open(fout,'w')

    GFF = GFFDB(gffdbfn)
    genome = Genome(seqfile)
    transcripts = GFF.features_of_type('mRNA')
    buffer = 4
    splice_halfwidth = readlen/2 + buffer
    transtable = string.maketrans("ATCG", "TAGC")
    for t in transcripts:
        chrom = t.chr
        strand = t.strand
        exons = list(GFF.children(t.id))
        if len(exons) == 1:
            continue # no splices, so no splice junctions

        if strand == '+':
            for i in range(len(exons)-1):
                splicestart1 = exons[i].stop - splice_halfwidth
                splicestop1 = exons[i].stop
                seq1 = genome.sequence(chrom,splicestart1,splicestop1)
                splicestart2 = exons[i+1].start
                splicestop2 = exons[i+1].start + splice_halfwidth + buffer
                seq2 = genome.sequence(chrom,splicestart2, splicestop2)
                description = '>%s%s%s%s%s' % (t.id, separator, exons[i].id, separator, exons[i+1].id)
                seq = 'NN' + seq1 + seq2 + 'NN'
                fout.write(description+'\n')
                fout.write(seq+'\n')  

        if strand == '-':
            for i in range(1,len(exons)):
                splicestart1 = exons[-i-1].stop - splice_halfwidth
                splicestop1 = exons[-i-1].stop
                seq1 = genome.sequence(chrom,splicestart1,splicestop1)
                splicestart2 = exons[-i].start
                splicestop2 = exons[-i].start + splice_halfwidth + buffer
                seq2 = genome.sequence(chrom,splicestart2, splicestop2)
                seq = 'NN' + seq1 + seq2 + 'NN'
                seq.translate(transtable)[::-1]
                description = '>%s%s%s%s%s' % (t.id, separator, exons[-i-1].id, separator, exons[-i].id)
                seq = 'NN' + seq1 + seq2 + 'NN'
                fout.write(description+'\n')
                fout.write(seq+'\n')  
       
        if close_fout:
            fout.close()

def merge_features(features, ignore_strand=False):
    """
    *features* is an iterable of features that you want to merge together.

    Will be converted into a list because it needs to sort the features.

    Returns an iterator of the merged features.

    *features* must be all on the same strand.  The new *featuretype* will be a 
    joining of the various featuretypes that were provided.

    If *ignore_strand* is True, strand will be forced to all '+' and no checking
    will be done on the input
    """
    # If it quacks like an iterator...then turn it into a list.
    if hasattr(features,'next'):
        features = list(features)
    
    # Quick function to sort by start position
    def start_pos(x):
        return x.start
    
    # Sort 'em by start position
    features.sort(key=start_pos)
    
    # Not sure how to deal with multiple strands...
    if not ignore_strand:
        strands = [i.strand for i in features]
        if len(set(strands))!= 1:
            raise NotImplementedError, 'Merging multiple strands not implemented yet'
        strand = strands[0]
    else:
        strand = '+'
    
    # ...or multiple chromosomes
    chroms = [i.chr for i in features]
    if len(set(chroms)) != 1:
        raise NotImplementedError, 'Merging multiple chromosomes not implemented yet'
    chrom = chroms[0]

    # If they were all exons, merged objects will have featuretype
    # 'merged_exon'.  If features included both exon and CDS, will have
    # featuretype 'merged_exon_CDS'.

    featuretypes = list(set([i.featuretype for i in features]))
    featuretypes.sort()
    featuretypes = '_'.join(featuretypes)
    featuretype = 'merged_%s' % featuretypes

    for i,feature in enumerate(features):
        if i == 0:
            current_merged_start = feature.start
            current_merged_stop = feature.stop
            continue
        # Does this feature start within the currently merged feature?...
        if feature.start <= current_merged_stop+1:
            # ...It starts within, so leave current_merged_start where it
            # is.  Does it extend any farther?
            if feature.stop >= current_merged_stop:
                # Extends further, so set a new stop position
                current_merged_stop = feature.stop
            else: 
                # If feature.stop < current_merged_stop, it's completely
                # within the previous feature.  Nothing more to do.
                continue
        else:
            # The start position is within the merged feature, so we're done with 
            # the current merged feature.  Prepare for output...

            merged_feature = GFFFeature(chr=feature.chr,
                                     source=None,
                                     featuretype=featuretype,
                                     start=current_merged_start,
                                     stop=current_merged_stop,
                                     value=None,
                                     strand=strand,
                                     phase=None,
                                     attributes=None)
            yield merged_feature

            # and we start a new one, initializing with this feature's start and stop.
            current_merged_start = feature.start
            current_merged_stop = feature.stop

    # need to yield the last one.
    #Feature = feature.__class__
    merged_feature = GFFFeature(chr=feature.chr,
                             source=None,
                             featuretype=featuretype,
                             start=current_merged_start,
                             stop=current_merged_stop,
                             value=None,
                             strand=strand,
                             phase=None,
                             attributes=None)
    yield merged_feature            

def intersect(features1, features2,ignore_strand=False):
    """
    Performs a pairwise intersection between 2 sets of features.
    """
    intersections = []
    if hasattr(features1,'next'):
        features2 = list(features1)
    if hasattr(features2,'next'):
        features2 = list(features2)

    for feature1 in features1:
        for feature2 in features2:
            if feature1.chr != feature2.chr:
                continue
            if not ((feature1.stop > feature2.start) and (feature1.start < feature2.stop)):
                continue

            maxstart = max( [feature1.start, feature2.start] )
            minstop = min( [feature1.stop, feature2.stop] )
            newfeature = GFFFeature(chr=feature1.chr,
                                     source=None,
                                     featuretype='intersection',
                                     start=maxstart,
                                     stop=minstop,
                                     value=None,
                                     strand=feature1.strand,
                                     phase=None,
                                     attributes=None)

            intersections.append( newfeature )
    for i in merge_features(intersections,ignore_strand=ignore_strand):
        i.featuretype = 'intersection'
        yield i

def create_gffdb(gfffn, dbfn):
    """
    Reads in a GFF3 file and constructs a database for use with downstream
    analysis.  See the GFFDB class in particular for an interface to this
    database. 

    This takes 2-3 min on a 100 MB GFF file for D. melanogaster.  This is a
    one-shot time investment, since once it's created using the database is
    quite fast.

    Note that the % complete for the first step is a percentage of the lines in
    the file.  Some GFF files have the entire sequence appended to the end
    which may cause the percentages to appear low.
    """

    # Calculate lines so you can display the percent complete
    if os.path.splitext(gfffn)[-1] == '.gz':
        f = gzip.open(gfffn)
    else:
        f = open(gfffn)
    nlines = 0.0
    for line in f:
        nlines += 1
    f.close()

    # This dictionary will contain the featurecounts of unlabeled features
    # for assigning unique IDs to features with no ID in their attributes.
    # Values will increment when a feature has no ID.
    #
    #   e.g., 
    #       featurecounts {'gene':42,
    #                      'three_prime_UTR':3,
    #                      'exon':408}
    featurecounts = {}


    # Create tables and indexes.
    conn = sqlite3.connect(dbfn)

    # Need this, otherwise you'll get annoying Unicode errors when reading in
    # GFF files that have special characters in their attributes.
    conn.text_factory = sqlite3.OptimizedUnicode

    c = conn.cursor()
    
    # Create the features table
    c.executescript('''
    CREATE TABLE features (
                            id text, 
                            chrom text, 
                            start int, 
                            stop int, 
                            strand text,
                            featuretype text,
                            value float, 
                            source text,
                            phase text,
                            attributes text,
                            primary key (id)
                          );
    CREATE TABLE relations (parent text, child text, level int, primary key(parent,child,level) );
    ''')

    # Parse the GFF file, populating the features tables and the relations table.
    t0 = time.time()
    counter = 0
    last_perc = 0
    for feature in GFFFile(gfffn,strvals=True):
        counter += 1
        perc_done = int(counter/nlines*100)
        if perc_done != last_perc:
            print '\rpopulating database with features and first-order parents: %s%%'%perc_done,
            sys.stdout.flush()
        last_perc = perc_done

        # If no ID, then assign one.
        if feature.id is None:
            featurecounts.setdefault(feature.featuretype,0)
            featurecounts[feature.featuretype] += 1
            new_id = 'unnamed_%s_%s' % (feature.featuretype, featurecounts[feature.featuretype]) 
            feature.add_attribute('ID',new_id)

        # actual insertion here
        c.execute('''
                  INSERT INTO features VALUES (?,?,?,?,?,?,?,?,?,?)
                  ''',(feature.id, 
                   feature.chr,
                   feature.start, 
                   feature.stop,
                   feature.strand,
                   feature.featuretype,
                   feature.value,
                   feature.source,
                   feature.phase,
                   feature._strattributes))

        try:
            parents = feature.attributes.Parent
            child = feature.id
        except AttributeError: 
            continue

        # add the first-level parent relation
        for parent in parents:
            c.execute('INSERT INTO relations VALUES (?,?,?)', (parent, child, 1))
    
    # show how much time it took
    t1 = time.time()
    print '(%ds)' % (t1-t0)
    t0 = time.time()

    # Index the database to speed up the second-level queries.
    print 'creating indexes: 0%',
    sys.stdout.flush()
    c.execute('CREATE INDEX ids ON features (id)')
    print '\rcreating indexes: 33%',
    sys.stdout.flush()
    c.execute('create index parentindex on relations (parent)')
    print '\rcreating indexes: 66%',
    sys.stdout.flush()
    c.execute('create index childindex on relations (child)')
    print '\rcreating indexes: 100%',
    sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    # OK, now go through everything again to find the second-order children
    t0 = time.time()
    c.execute('select id from features')

    # create 2 more cursors so you can iterate over one while querying on the other's iteration
    c2 = conn.cursor()
    c3 = conn.cursor()

    # For each feature in the features table, query the the relations table
    counter = 0
    grandchild_count = 0
    last_perc = 0
    tmp = tempfile.mktemp()
    fout = open(tmp,'w')
    for parent in c:
        parent = parent[0] # first thing in the list is the ID
        counter += 1
        perc_done = int(counter / nlines * 100)
        if perc_done != last_perc:
            print '\rquerying database for grandchildren: %s%%'%perc_done,
            sys.stdout.flush()
        last_perc = perc_done    

        # Here we get the first-level child from the initial import.   This
        # data is contained in the "Parent=" attribute of each GFF feature.
        c2.execute('select child from relations where parent = ? and level=1',(parent,))
        for child in c2:
            child = child[0]
            c3.execute('select child from relations where parent = ? and level=1',(child,))
            for grandchild in c3:
                grandchild_count += 1
                grandchild = grandchild[0]
                fout.write('%s\t%s\n' % (parent,grandchild))
    fout.close()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    t0 = time.time()
    counter = 0
    last_perc = 0
    grandchild_count = float(grandchild_count)
    for line in open(tmp):
        counter += 1
        perc_done = int(counter / grandchild_count * 100)
        if perc_done != last_perc:
            print '\rimporting grandchildren: %s%%' % perc_done,
            sys.stdout.flush()
        last_perc = perc_done
        parent,child = line.strip().split('\t')
        c.execute('insert or ignore into relations values (?,?,?)', (parent, child, 2))
    t1 = time.time()
    print '(%ds)' % (t1-t0)
    t0 = time.time()
    print 're-creating indexes',
    c.execute('drop index childindex')
    c.execute('drop index parentindex')
    c.execute('create index parentindex on relations (parent)')
    print '\rre-creating indexes: 50%',
    sys.stdout.flush()
    c.execute('create index childindex on relations (child)')
    c.execute('create index starts on features(start)')
    c.execute('create index stops on features(stop)')
    c.execute('create index startstrand on features(start, strand)')
    c.execute('create index stopstrand on features(stop,strand)')
    c.execute('create index featuretypes on features(featuretype)')
    print '\rre-creating indexes: 100%',
    sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    os.remove(tmp)
        

def create_gtfdb(gtffn, dbfn):
    """
    Reads in a GTF file and constructs a database for use with downstream
    analysis.  The GTFDB class will work as an interface to this database. 

    Assume that for a CDS, the attribute gene_id specifies the "grandparent" gene
    and the attribute transcript_id specifies the "parent" transcript.

    This also implies that the gene_id is the parent of the transcript_id.

    Thus, the relations table is filled out first.
    """
    # Not all GTF files contain "exon_number" attributes.  We need these in
    # order to construct unique IDs for features. This dictionary will hold the
    # current highest count for each featuretype and for each gene.
    #
    # So feature_counts will look something like:
    #
    #     feature_counts = {'exon':{'gene1':1,
    #                               'gene2':8,},
    #                       '5UTR':{'gene1':1,
    #                               'gene2':1}
    #                      }

    feature_counts = {}


    # Calculate lines so you can display the percent complete
    f = open(gtffn)
    nlines = 0.0
    for line in f:
        nlines += 1
    f.close()

    # Create tables and indexes.
    conn = sqlite3.connect(dbfn)
    c = conn.cursor()
    c.executescript('''
    CREATE TABLE features (
                            id text primary key, 
                            chrom text, 
                            start int, 
                            stop int, 
                            strand text,
                            featuretype text,
                            value float, 
                            source text,
                            phase text,
                            attributes text
                          );
    CREATE TABLE relations (parent text, child text, level int, primary key (parent, child));
    ''')

    # Parse the GTF file, populating the features tables and the relations table.
    t0 = time.time()
    counter = 0
    last_perc = 0
    for feature in GTFFile(gtffn,strvals=True):
        counter += 1
        perc_done = int(counter/nlines*100)
        if perc_done != last_perc:
            print '\rpopulating database with features and first-order parents: %s%%'%perc_done,
            sys.stdout.flush()
        last_perc = perc_done

        parent = feature.attributes.transcript_id
        grandparent = feature.attributes.gene_id
        
        # Initialize counts to 1 and get the feature number which will become
        # the feature ID for the database.
        #
        # TODO: maybe this should be a uniquely-identifiable number (e.g.,
        # gene-start-stop?, rather than an integer that is dependent on the
        # location in the source file?
        feature_counts.setdefault(feature.featuretype,{}).setdefault(grandparent,0)
        feature_counts[feature.featuretype][grandparent] += 1
        feature_number = feature_counts[feature.featuretype][grandparent]
        
        # Note: with the new add_attribute() method adding to the string as
        # well, we don't want something that permanent . . . so just stick with
        # this just-for-the-database "ID" instead of
        # "feature.add_attribute('ID',ID)"
        ID = '%s:%s:%s-%s' % (feature.featuretype, parent, feature.start, feature.stop)
        

        # If it's an exon, its attributes include its parent transcript
        # and its 'grandparent' gene.  So we can insert these
        # relationships into the relations table now.

        # Note that the table schema has (parent,child) as a primary
        # key, so the INSERT OR IGNORE won't add multiple entries for a
        # single (parent,child) relationship

        # The gene has a grandchild exon
        c.execute('''
        INSERT OR IGNORE INTO relations VALUES (?,?,?)
        ''', (grandparent, ID, 2))

        # The transcript has a child exon
        c.execute('''
        INSERT OR IGNORE INTO relations VALUES (?,?,?)
        ''', (parent, ID, 1))

        # The gene has a child transcript
        c.execute('''
        INSERT OR IGNORE INTO relations VALUES (?,?,?)
        ''', (grandparent, parent, 1))

        # Insert the feature into the features table.
        c.execute('''
                  INSERT OR IGNORE INTO features VALUES (?,?,?,?,?,?,?,?,?,?)
                  ''',(ID, 
                   feature.chr,
                   feature.start, 
                   feature.stop,
                   feature.strand,
                   feature.featuretype,
                   feature.value,
                   feature.source,
                   feature.phase,
                   feature._strattributes))
    conn.commit()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    # Now that the entries in the GTF file have been added, we can index
    # the database to speed up the second-level queries we need to do to
    # identify genes.
    t0 = time.time()
    print 'creating indexes: 0%',
    sys.stdout.flush()
    c.execute('DROP INDEX IF EXISTS ids')
    c.execute('CREATE INDEX ids ON features (id)')
    print '\rcreating indexes: 33%',
    sys.stdout.flush()
    c.execute('DROP INDEX IF EXISTS parentindex')
    c.execute('create index parentindex on relations (parent)')
    print '\rcreating indexes: 66%',
    sys.stdout.flush()
    c.execute('drop index if exists childindex')
    c.execute('CREATE INDEX childindex ON relations (child)')
    print '\rcreating indexes: 100%',
    sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    
    # GTF files only describe exons, so genes and transcripts are described
    # implicitly.  With all the exons imported into the database, we can
    # now do queries to find the limits of genes and transcripts . . . 

    # Get an idea of how many parents we have to get through
    t0 = time.time()
    c.execute("SELECT COUNT(DISTINCT parent) FROM relations")
    nparents = float(c.fetchone()[0])
    
    # These results will include transcripts and genes.
    c.execute("""
              SELECT DISTINCT parent FROM relations
              """)
    tmp = tempfile.mktemp()
    fout = open(tmp,'w')
    counter = 0
    c2 = conn.cursor()
    for parent in c:

        # Some feedback...
        counter += 1
        perc_done = int(counter/nparents*100)
        if perc_done != last_perc:
            print '\rquerying database for start/stop positions of genes and transcripts: %s%%'%perc_done,
            sys.stdout.flush()
        last_perc = perc_done

        # Find the start and stop limits of this parent's children
        parent = parent[0]
        c2.execute("""
                   select min(start), max(stop), level, strand, chrom FROM features 
                   JOIN relations ON
                   features.ID = relations.child
                   WHERE
                   parent = ?
                   """, (parent,))
        
        # For testing to make sure you only get one level back.
        #child_limits = c2.fetchall()
        #assert len(child_limits) == 1
        #child_limits = child_limits[0]

        child_limits = c2.fetchone()
        start,end,level,strand, chrom = child_limits

        # The strategy here is to write parents to file, and later read the
        # file back into the database so that we don't keep writing to the
        # database while we're in the middle of consuming a query results
        # iterator...

        # Using some assumptions we can make some shortcuts to determine
        # what featuretype this parent really is.  Since the features table
        # only contains exons, and we're joining on the features table,
        # only exon children or only grandchildren will be returned (this
        # can be verified by the commented-out test code above).  If only
        # grandchildren were returned (i.e., level=2) then the parent is a
        # gene.

        # In addition, we need to create a fake attributes string so that
        # later, upon accessing the database, the GTFDB wrapper will know
        # how to assign an ID.  This is sort of a hack in order to maintain
        # cleaner class inheritance between GFFFeatures and GTFFeatures.

        # Since the relations table only has transcript children in it, any parents
        # at level 1 are mRNAs.  
        #
        # WARNING: this does NOT account for non-coding RNAs -- they will still
        # be called "mRNA"
        if level == 1:
            featuretype = 'mRNA'
            attributes = 'transcript_id "%s"; ' % parent

        # On the other hand, level 2 parents means that the parent is a gene.
        if level == 2:
            featuretype = 'gene'
            attributes = 'gene_id "%s"; ' % parent


        if level is None:
            print 'WARNING: got back nothing good from db for %s' % parent
            featuretype='None'

        fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (parent,chrom,
                                                     start,end,strand,
                                                     featuretype,attributes))
    fout.close()
    

    # Now that the file has been created and we're done querying for the
    # parents, slurp in the file and insert everything into the db.
    fin = open(fout.name)
    for line in fin:
        c.execute("""
                  INSERT OR IGNORE INTO features (id, chrom, start, stop,
                  strand, featuretype, attributes) VALUES (?,?,?,?,?,?,?)
                  """, line.strip().split('\t'))

    conn.commit()

    # clean up the tempfile.
    os.remove(tmp)

    # show how much time it took
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    # Re-create indexes since we just added all those parents. Prob only
    # need to re-create the first one though.
    t0 = time.time()
    print 're-creating indexes: 0%',
    sys.stdout.flush()
    c.execute('DROP INDEX ids')
    c.execute('CREATE INDEX ids ON features (id)')
    print '\rre-creating indexes: 33%',
    sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    print '(%ds)' % (t1-t0)

    conn.commit()
        
    del feature_counts

class GFFDB:
    featureclass = GFFFeature
    add_id = ''
    def __init__(self,db_fn):
        """*db_fn* is the filename of the database to use"""
        self.db_fn = db_fn
        self.conn = sqlite3.connect(db_fn)
        self.conn.text_factory = str

    def __getitem__(self,id):
        c = self.conn.cursor()
        c.execute('''
                  SELECT %s chrom, source, featuretype, start, stop, value, strand, phase, attributes from features where id = ?
                  ''' % self.__class__.add_id, (id,))
        results = c.fetchall()
        assert len(results) == 1, len(results)
        return self.__class__.featureclass(*results[0])

    def features_of_type(self, featuretype, chrom=None, start=None, stop=None, strand=None):
        """
        Returns an iterator of GFFFeature objects that have the feature type
        *featuretype*.  You can optionally specify chrom, strand, start, and
        stop filters.

        For example, to get all exons within a range::

            self.features_of_type('exon', chrom='chr3L', start=11100, stop=12200)

        Any of the filters can be empty::

            # exons on any chromosome after the first 1kb
            self.features_of_type('exon',start=1000)

            # genes on + strand of chrX
            self.features_of_type('gene',chrom='chrX', strand='+')

        For more complicated stuff, (e.g., filtering by multiple chromosomes)
        you'll probably want to use self.execute() directly and write your own
        queries.
        """
        filter_clause = ''
        if chrom is not None:
            filter_clause += ' AND chrom = "%s"' % chrom

        if start is not None:
            filter_clause += ' AND start >= %s' % start
        if stop is not None:
            filter_clause += ' AND stop <= %s' % stop
        if strand is not None:
            filter_clause += ' AND strand = "%s"' % strand

        c = self.conn.cursor()
        c.execute('''
                  SELECT 
                  id,chrom, source, featuretype, start, stop, value, strand,
                  phase, attributes
                  FROM 
                  features 
                  WHERE featuretype = ?
                  %s
                  ORDER BY start
                  ''' % filter_clause , (featuretype,))
        for i in c:
            yield self[i[0]]
    
    def features(self):
        """
        Returns an iterator of the different feature types found in the database.
        """
        c = self.conn.cursor()
        c.execute('''
                  SELECT DISTINCT featuretype from features
                  ''')
        for i in c:
            yield i[0]

    def all(self):
        """
        Returns an iterator of ALL features in the database.
        """
        c = self.conn.cursor()
        c.execute('''
                  SELECT
                  %s chrom, source, featuretype, start, stop, value, strand,
                  phase, attributes
                  FROM
                  features 
                  ''' % self.__class__.add_id)
        for i in c:
            yield self.__class__.featureclass(*i)

    def closest_feature(self, chrom, pos, featuretype='gene', strand=None, ignore=None, direction=None):
        """
        Returns the distance and ID of the closest TSS to the coordinate. Strand optional.
        """
        
        # e.g., AND id != FBgn0001 AND id != FBgn0002
        ignore_clause = ''
        if ignore is not None:
            if type(ignore) is str:
                ignore = [ignore]
            for i in ignore:
                ignore_clause += ' AND id != "%s" ' % i
        if (strand is None) and (direction is not None):
            raise ValueError, 'Strand must be specified if direction is specified'

        # Mini-examples above each query show the position as an "x" and
        # promoters to indicate genes on (+) or (-) strands

        #  ---->                ----->
        #  |               X    |
        #  ^ finds this one
        if strand == '+' and direction == 'upstream':
            c = self.conn.cursor()
            c.execute('''
            SELECT %s-start as AAA,id FROM features 
            WHERE featuretype = ?
            AND chrom = ?
            AND strand = "+"
            AND start < ?
            %s
            ORDER BY AAA LIMIT 1
            ''' % (pos, ignore_clause), (featuretype,chrom,pos))
            closest = c.fetchone()
            return closest
        
        #  ---->                 ----->
        #  |       X             |
        #                        ^ finds this one
        if strand == '+' and direction == 'downstream':
            c = self.conn.cursor()
            c.execute('''
            SELECT start-%s as AAA,id FROM features 
            WHERE featuretype = ?
            AND chrom = ?
            AND strand = "+"
            AND start > ?
            %s
            ORDER BY AAA LIMIT 1
            ''' % (pos, ignore_clause), (featuretype,chrom,pos))
            closest = c.fetchone()
            return closest

        #  <----                   <----
        #      |    X                  |
        #                              ^ finds this one
        if strand == '-' and direction == 'upstream':
            c = self.conn.cursor()
            c.execute('''
            SELECT stop-%s as AAA,id FROM features 
            WHERE featuretype = ?
            AND chrom = ?
            AND strand = "-"
            AND stop > ?
            %s
            ORDER BY AAA LIMIT 1
            ''' % (pos, ignore_clause), (featuretype,chrom,pos))
            closest = c.fetchone()
            return closest

        #  <----               <----
        #      |           X       |
        #      ^ finds this one
        if strand == '-' and direction == 'downstream':
            c = self.conn.cursor()
            c.execute('''
            SELECT %s-stop as AAA,id FROM features 
            WHERE featuretype = ?
            AND chrom = ?
            AND strand = "-"
            AND stop < ?
            %s
            ORDER BY AAA LIMIT 1
            ''' % (pos, ignore_clause), (featuretype,chrom,pos))
            closest = c.fetchone()
            return closest

        if direction is None:
            #  ---->     <--         ---->
            #  |           |   X     |
            #                        ^ finds this one
            if strand == '+':
                c = self.conn.cursor()
                c.execute('''
                SELECT abs(start-%s) as AAA,id FROM features 
                WHERE featuretype = ?
                AND chrom = ?
                AND strand = "+"
                %s
                ORDER BY AAA LIMIT 1
                ''' % (pos, ignore_clause), (featuretype,chrom))
                closest = c.fetchone()
                return closest

            #  <----     --->  <----
            #      |     |   X     |
            #                      ^ finds this one
            if strand == '-':
                c = self.conn.cursor()
                c.execute('''
                SELECT abs(stop-%s) as AAA,id FROM features
                WHERE featuretype = ?
                AND chrom = ?
                AND strand = "-"
                %s
                ORDER BY AAA LIMIT 1
                ''' % (pos, ignore_clause), (featuretype,chrom))
                closest = c.fetchone() 
                return closest


            # If strand not specified, then find the closest one on either strand
            #  <----             -->      <----
            #      |         X   |            |
            #                    ^ finds this one
            if strand is None:
                c = self.conn.cursor()
                c.execute('''
                SELECT abs(start-%s) as AAA,id FROM features 
                WHERE featuretype = ?
                AND chrom = ?
                AND strand = "+"
                %s
                ORDER BY AAA LIMIT 1
                ''' % (pos, ignore_clause), (featuretype,chrom))
                closest_plus = c.fetchone()

                c.execute('''
                SELECT abs(stop-%s) as AAA,id FROM features
                WHERE featuretype = ?
                AND chrom = ?
                AND strand = "-"
                %s
                ORDER BY AAA LIMIT 1
                ''' % (pos, ignore_clause), (featuretype,chrom))
                closest_minus = c.fetchone() 
                both = [closest_minus,closest_plus]
                both.sort()
                return both[0]

    def overlapping_features(self,chrom,start,stop,featuretype=None,strand=None, completely_within=False):
        """
        Returns an iterator of features of type *featuretype* that overlap the
        feature with ID *id*. If *featuretype* is None (default), then all
        features will be returned.

        If *ignore_strand* is True, features of any strand will be returned.
        If *ignore_strand* is False (default) then only return overlapping
        features on the same strand as *id*.

        If *completely_within* is True, a feature needs to be completely within
        the boundaries of *id*.  If False (default), features that extend out
        on either side will be included as well.
        """
        
        if strand is None:
            strand_clause = ''
        else:
            strand_clause = ' AND strand = "%s"' % strand
        
        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND featuretype = "%s"' % featuretype

        if completely_within:
            within_clause = ' AND ((start BETWEEN %s AND %s) AND (stop BETWEEN %s AND %s))' % (start,stop, start,stop)
        else:
            within_clause = ' AND ((start BETWEEN %s AND %s) OR (stop BETWEEN %s AND %s))' % (start,stop, start,stop)
        c = self.conn.cursor()
        c.execute('''
        SELECT %s chrom,source,featuretype,start,stop,value,strand,phase,attributes
        FROM features WHERE
        chrom = ?
        %s
        %s
        %s
        ''' % (self.add_id, strand_clause, featuretype_clause, within_clause), (chrom,))
        for i in c:
            yield self.__class__.featureclass(*i)

    def merge_features(self, features, ignore_strand=False):
        """
        *features* is an iterable of features that you want to merge together.

        Will be converted into a list because it needs to sort the features.

        Returns an iterator of the merged features.

        *features* must be all on the same strand.  The new *featuretype* will be a 
        joining of the various featuretypes that were provided.

        If *ignore_strand* is True, strand will be forced to all '+' and no checking
        will be done on the input
        """
        # If it quacks like an iterator...then turn it into a list.
        if hasattr(features,'next'):
            features = list(features)
        
        # Quick function to sort by start position
        def start_pos(x):
            return x.start
        
        # Sort 'em by start position
        features.sort(key=start_pos)
        
        # Not sure how to deal with multiple strands...
        if not ignore_strand:
            strands = [i.strand for i in features]
            if len(set(strands))!= 1:
                raise ValueError, 'Specify ignore_strand=True to force merging of multiple strands'
            strand = strands[0]
        else:
            strand = '+'
        
        # ...or multiple chromosomes
        chroms = [i.chr for i in features]
        if len(set(chroms)) != 1:
            raise NotImplementedError, 'Merging multiple chromosomes not implemented yet'
        chrom = chroms[0]

        # If they were all exons, merged objects will have featuretype
        # 'merged_exon'.  If features included both exon and CDS, will have
        # featuretype 'merged_exon_CDS'.

        featuretypes = list(set([i.featuretype for i in features]))
        featuretypes.sort()
        featuretypes = '_'.join(featuretypes)
        featuretype = 'merged_%s' % featuretypes

        for i,feature in enumerate(features):
            if i == 0:
                current_merged_start = feature.start
                current_merged_stop = feature.stop
                continue
            # Does this feature start within the currently merged feature?...
            if feature.start <= current_merged_stop+1:
                # ...It starts within, so leave current_merged_start where it
                # is.  Does it extend any farther?
                if feature.stop >= current_merged_stop:
                    # Extends further, so set a new stop position
                    current_merged_stop = feature.stop
                else: 
                    # If feature.stop < current_merged_stop, it's completely
                    # within the previous feature.  Nothing more to do.
                    continue
            else:
                # The start position is within the merged feature, so we're done with 
                # the current merged feature.  Prepare for output...
                Feature = feature.__class__
                merged_feature = Feature(chr=feature.chr,
                                         source=None,
                                         featuretype=featuretype,
                                         start=current_merged_start,
                                         stop=current_merged_stop,
                                         value=None,
                                         strand=strand,
                                         phase=None,
                                         attributes=None)
                yield merged_feature

                # and we start a new one, initializing with this feature's start and stop.
                current_merged_start = feature.start
                current_merged_stop = feature.stop

        # need to yield the last one.
        Feature = feature.__class__
        merged_feature = Feature(chr=feature.chr,
                                 source=None,
                                 featuretype=featuretype,
                                 start=current_merged_start,
                                 stop=current_merged_stop,
                                 value=None,
                                 strand=strand,
                                 phase=None,
                                 attributes=None)
        yield merged_feature            

    def interfeatures(self,features):
        """
        Given an iterable of *features*, returns an iterator of new features
        defined by the intervening space between each consecutive feature and
        the one before it.

        If you provide N features, you'll get back N-1 intervening features.
        
        This is a purposefully naive method that does NOT do sorting or merging
        -- it simply returns the intervening space between the features provided.

        So if *features* contains the exons of a transcript, this method will
        return the introns.

        If the features returned are not sorted, you may get overlapping
        results.

        Example for getting all the 'non-exonic' space::
            
            # merge_features needs a single chrom and single strand
            exons = G.features_of_type('exon',chrom='chrX',strand='+')
            merged_exons = G.merge_features(exons)
            non_exonic = G.interfeatures(merged_exons)
        """
        for i,feature in enumerate(features):
            if i == 0:
                interfeature_start = feature.stop
                last_feature = feature
                continue
            interfeature_stop = feature.start
            featuretype = 'inter_%s_%s' % (last_feature.featuretype, feature.featuretype)
            assert last_feature.strand == feature.strand
            assert last_feature.chr == feature.chr
            strand = last_feature.strand
            chr = last_feature.chr

            # Shrink
            interfeature_start += 1
            interfeature_stop -= 1

            yield self.featureclass(chr=chr,
                                    source=None,
                                    featuretype=featuretype,
                                    start=interfeature_start,
                                    stop=interfeature_stop,
                                    value=None,
                                    strand=strand,
                                    phase=None,
                                    attributes=None)
            interfeature_start = feature.stop

    def chromosomes(self):
        """
        Returns a list of chromosomes in the database
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT DISTINCT chrom FROM features
        ''')
        return [i[0] for i in c.fetchall()]

    def strands(self):
        """
        Returns a list of strands in the database
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT DISTINCT strand FROM features
        ''')
        return [i[0] for i in c.fetchall()]
 
    def attach_children(self,feature):
        '''
        Adds a new attribute to GFFFeature objects, children.  This is a list of GFFFeatures
        that are the children of the transcript.  This also adds "introns" as children.

        Returns an iterator, of length 1 if *id* was a transcript, or of length equal to number of
        child transcripts'''
         
        def sortfunc(x):
            return x.start
        
        if feature.featuretype == 'gene':
            raise NotImplementedError,'not sure how to deal with genes yet'
        
        children = list(self.children(feature.id))
        
        # CDSs are duplicated as 'exons'.  So remove the exons from the
        # list of children.
        
        # annoying: some, but not all, exons are duplicated as CDSs.  Remove
        # the exons that have the same start/stop as CDSs.
             
        CDScoords = [(i.start,i.stop) for i in children if i.featuretype=='CDS']
        new_children = []
        for child in children:
            if child.featuretype == 'exon':
                if (child.start,child.stop) in CDScoords:
                    continue
            new_children.append(child)

        children = new_children

        if feature.featuretype == 'mRNA':
            children_for_intron_calc = [i for i in children if i.featuretype != 'CDS']
        else:
            children_for_intron_calc = children

        children_for_intron_calc.sort(key=sortfunc)
        # construct introns:
        exons = [ (i.start,i.stop) for i in children_for_intron_calc ]
        for i in range(len(exons)-1):
            thisstart,thisstop = exons[i]
            nextstart,nextstop = exons[i+1]
            intron = GFFFeature(chr=feature.chr,
                                source='imputed',
                                featuretype='intron',
                                start=thisstop+1,
                                stop=nextstart-1,
                                value=None,
                                strand=feature.strand,
                                phase=None,
                                attributes=None)
            children.append( intron )
            
        # sort children by increasing start position.
        children.sort(key=sortfunc)
        feature.children = children
        return feature
   
    def children(self,id,level=1,featuretype=None):
        '''Returns an iterator of the children *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified.'''
        cursor = self.conn.cursor()

        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND features.featuretype = "%s"' % featuretype

        cursor.execute('''
            SELECT DISTINCT
            %s chrom, source, featuretype, start, stop, value, strand, phase, attributes 
            FROM features JOIN relations 
            ON relations.child = features.id
            WHERE relations.parent = ? 
            AND relations.level = ? 
            %s
            ORDER BY start''' % (self.__class__.add_id, featuretype_clause),(id,level))
        for i in cursor:
            yield self.__class__.featureclass(*i)

    def parents(self,id,level=1,featuretype=None):
        '''Returns an iterator of the parents *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified.'''
        cursor = self.conn.cursor()
        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND features.featuretype = "%s"' % featuretype
        cursor.execute('''
            SELECT DISTINCT
            %s chrom, source, featuretype, start, stop, value, strand, phase, attributes 
            FROM features JOIN relations 
            ON relations.parent = features.id
            WHERE relations.child = ? 
            AND relations.level = ? 
            %s
            ORDER BY start''' % (self.__class__.add_id, featuretype_clause),(id,level))
        for i in cursor:
            yield self.__class__.featureclass(*i)

    def execute(self, query):
        """
        Returns an iterator of results from any arbitrary query passed in.

        No type conversion is done, so you won't get GFFFeature objects back,
        just tuples of strings.
        """
        cursor = self.conn.cursor()
        cursor.execute(query)
        for i in cursor:
            yield i

    def refFlat(self, gene):
        """Writes the gene out as a RefFlat format:
            
            geneName altname chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds
        
        Caveats: 

            cdsStart and cdsEnd are the min and max positions, respectively, of all CDSs in the gene.
            So a particular isoform's CDS may not actually start on the gene-wide minimum CDS position.

            Assumes that there was an attribute in the GFF file called 'Name'.  This will then be
            found in the gene.attributes.Name attribute.
        """
        
        geneID = gene.id
        txStart = gene.start
        txEnd = gene.stop
        exons = []
        cdss = []
        exon_count = 0
        for i in self.children(geneID,2):
            if i.featuretype == 'CDS':
                cdss.append(i)
            if i.featuretype == 'exon':
                exons.append(i)
        
        if len(cdss) == 0:
            cdsStart = 'NA'
            cdsEnd = 'NA'
        else:
            cdsStart = min([i.start for i in cdss])
            cdsEnd = max([i.stop for i in cdss])

        def exon_sort(e):
            return e.start

        exons.sort(key=exon_sort)


        exonStarts = ','.join([str(i.start) for i in exons])
        exonEnds = ','.join([str(i.stop) for i in exons])

        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene.attributes.Name[0], 
                                                                 gene.id,
                                                                 gene.chr, 
                                                                 gene.strand, 
                                                                 gene.start, 
                                                                 gene.stop,
                                                                 cdsStart,
                                                                 cdsEnd,
                                                                 len(exons),
                                                                 exonStarts,
                                                                 exonEnds)
        return line

    def count_features_of_type(self,featuretype):
        """
        More memory efficient than the functionally equivalent
        ::
            len(list(self.features_of_type(featuretype)))
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT count() FROM features
        WHERE featuretype = ?''',(featuretype,))
        results = c.fetchone()
        if results is not None:
            results = results[0]
        return results
    
    def promoter(self, id, dist=1000, truncate_at_next_feature=None, bidirectional=True):
        """
        Returns a new GFFFeature of featuretype "promoter", with the definition
        of the promoter guided by the kwargs described below.
        
        *dist* (default 1000) is the distance in bp from TSS that you want to
        include.  TSS is considered the feature start position if the feature
        is on the plus strand, or the feature stop position if on the minus
        strand. The application of *dist* to getting a promoter can be changed
        by the other kwargs below:

        If *bidirectional* is True (default), then include *dist* bp on BOTH
        sides of the TSS.  Otherwise, only include the region *dist* bp 5' of
        the TSS.
        
        If *truncate_at_next_feature* (None by default) is a featuretype (e.g.,
        'gene' or 'mRNA') then first try and go up to *dist* bp upstream from
        feature start.  If there is another feature of the specified
        featuretype within *dist* upstream on either strand, then the promoter
        region returned will be truncated to stop at this feature.
        
        If *truncate_at_next_feature* is not None AND *bidirectional* is True,
        then the downstream distance will not be truncated.

        Note that without supplying the chromosome limits, it's impossible to
        know whether the returned promoter runs off the end of the chromosome.
        The lower boundary is accounted for (i.e., no negative coordinates) but
        checking the upper boundary is the responsibility of the calling
        function.
        """
        c = self.conn.cursor()
        feature = self[id]
        
        if truncate_at_next_feature is not None:
            # Locate closest feature.  Different strands have different queries:
            #
            # for (+) strand, look for other features' stops that are less than
            # to this features's start and find out how far away it is:
            if feature.strand == '+':
                c.execute('''
                SELECT min(%s-stop)
                FROM features WHERE
                featuretype = ?
                AND chrom = ?
                AND id != ?
                AND stop < ?
                ''' % (feature.start),(truncate_at_next_feature, 
                                       feature.chr, 
                                       feature.id, 
                                       feature.start))
                closest_feature_dist = c.fetchone()[0]
                if (closest_feature_dist < dist) and (closest_feature_dist is not None):
                    dist = closest_feature_dist
            
            # for (-) strand, look for other features' starts that are greater
            # than this feature's stop and find out how far away it is:
            if feature.strand == '-':
                c.execute('''
                SELECT min(start-%s)
                FROM features WHERE
                featuretype = ?
                AND chrom = ?
                AND id != ?
                AND start > ?
                ''' % (feature.stop),(truncate_at_next_feature, 
                                      feature.chr, 
                                      feature.id, 
                                      feature.stop))
                closest_feature_dist = c.fetchone()[0]
                if (closest_feature_dist < dist) and (closest_feature_dist is not None):
                    dist = closest_feature_dist

        if feature.strand == '+':
            TSS = feature.start
            upstream = TSS - dist
            downstream = TSS + dist

        if feature.strand == '-':
            TSS = feature.stop
            upstream = TSS + dist
            downstream = TSS - dist
        
        # Negative coords not allowed; truncate to beginning of chrom (note the
        # incrementing of start and decrementing of stop below, so this will
        # become 1 instead of 0)
        if upstream < 0:
            upstream = 0

        if downstream < 0:
            downstream = 0

        if bidirectional:
            coords = [upstream,downstream]
        else:
            coords = [upstream, TSS]

        coords.sort()
        start,stop = coords

        # shrink by one on each side so there's no overlap
        start = start + 1
        stop = stop -1 

        # Create a new GFFFeature object (or GTFFeature) to return
        promoter = self.__class__.featureclass(chr=feature.chr,
                                              source='imputed',
                                              featuretype='promoter',
                                              start=start,
                                              stop=stop,
                                              strand=feature.strand,
                                              value=None,
                                              phase=None,
                                              attributes=None)
        promoter.add_attribute('ID','promoter:'+feature.id)
        return promoter

    def random_feature(self,featuretype=None):
        """
        Chooses a random feature from the database.  Useful for testing or
        experimenting with the module.  Specify a *featuretype* to restrict results
        to that feature type.
        
        Idea from here:
            
            http://www.mail-archive.com/sqlite-users@sqlite.org/msg14657.html
        """

        c = self.conn.cursor()
        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = 'featuretype = "%s" AND ' % featuretype
        c.execute('''
        SELECT %s chrom, source, featuretype, start, stop, value, strand, phase, attributes 
        FROM features
        WHERE 
        %s
        rowid >= abs(random()) %% (SELECT MAX(rowid) FROM features) LIMIT 1
        ''' % (self.__class__.add_id, featuretype_clause))
        results = c.fetchone()
        return self.__class__.featureclass(*results)
        

    def coding_genes(self):
        """
        Returns an iterator of genes that contain CDSs. Useful for if you want
        to exclude tRNA, various ncRNAs, etc, since they are also annotated
        with featuretype "gene".
        """
        for g in self.features_of_type('gene'):
            for grandchild in self.children(g.id, level=2):
                if grandchild.featuretype == 'CDS':
                    yield g
                    break
    
    def n_gene_isoforms(self, geneID):
        """
        Returns the number of isoforms that this gene has.
        """
        n = 0
        for i in self.children(geneID, level=1, featuretype='mRNA'):
            n += 1
        return n

    def n_exon_isoforms(self, exonID):
        """
        Returns the number of isoforms that this exon is found in.
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT count() FROM relations 
        JOIN features 
        ON relations.parent=features.id
        WHERE relations.child=?
        AND
        relations.level=1
        AND 
        features.featuretype="mRNA"
        ''', (exonID,))
        return c.fetchone()[0]

    def exons_gene(self,exonID):
        """
        Returns the ID of the exon's parent gene.  Fast, single-purpose
        method that doesn't do the type conversion or sorting of
        self.parents().
        """
        c = self.conn.cursor()
        c.execute('''
        SELECT parent FROM relations
        JOIN features 
        ON relations.parent=features.id
        WHERE child=?
        AND 
        relations.level=2
        AND
        features.featuretype="gene"
        ''', (exonID,))
        return c.fetchone()[0]


class GTFDB(GFFDB):
    featureclass = GTFFeature
    add_id = 'id,'

    def three_prime_UTR(self,id):
        """
        Returns the 3' UTR of a transcript (if *id*.featuretype == 'mRNA')
        or of all children transcripts (if *id*.featuretype == 'gene').

        Calculates the 3'UTR by doing coordinate differences between child
        exon and CDS features.
        """

        feature = self[id]
        if feature.featuretype == 'gene':
            children = self.children(id,level=2)
        if feature.featuretype == 'mRNA':
            children = self.children(id,level=1)
        elif feature.featuretype not in ['gene','mRNA']:
            raise NotImplementedError, "3' UTR not implemented for features of type %s" % feature.featuretype

        CDSs = []
        exons = []
        for child in children:
            if child.featuretype == 'CDS':
                CDSs.append(child)
            elif child.featuretype == 'exon':
                exons.append(child)
        
        # Now we have sorted features.
        UTRs = []
        for exon in exons:
            for CDS in CDSs:
                if exon.start <= CDS.start <= exon.stop:
                    if (exon.start == CDS.start) and (exon.stop == CDS.stop):
                        # exactly matches up, so not a UTR.
                        continue
                    else:
                        # partial overlap, so have to get the limits.
                        if exon.stop == CDS.stop:
                            UTRstart = exon.start
                            UTRstop = CDS.stop
                        elif exon.start == CDS.start:
                            UTRstart = CDS.stop
                            UTRstop = exon.stop
                else:
                    # If no overlap at all with a CDS, then the whole thing
                    # is a UTR.
                    UTRstart = exon.start
                    UTRstop = exon.stop
                
                UTR = copy.deepcopy(exon)
                UTR.id = UTR.id.replace('exon','UTR')
                UTR.start = UTRstart
                UTR.stop = UTRstop
                UTR.featuretype = 'UTR'
                UTRs.append(UTR)
        for i in UTRs:
            yield i

if __name__ == "__main__":
    G = GTFDB('/home/ryan/data/solexa/meta-solexa/dm3-chr.GTF.db')
