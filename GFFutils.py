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
import logging

logging.basicConfig(level=logging.DEBUG)

class FeatureNotFoundError(Exception):
    """
    Error to be raised when an ID is not in the database.
    """
    def __init__(self, feature_id):
        Exception.__init__(self)
        self.feature_id = feature_id
    def __str__(self):
        return self.feature_id

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
            self._attrs = []  # will hold a list of attributes added to the
                              # object.

    def __init__(self, chrom=None, source=None, featuretype=None, start=None, stop=None,
                 value=None, strand=None, phase=None, attributes=None, name=None, strvals=False):
        """
        Represents a line in a GFF file.

        *chrom*       : Chromosome
        *source*      : Source of the data
        *featuretype* : type of the feature ('gene', 'mRNA', 'CDS', etc)
        *start*       : Start position on the chromosome
        *stop*        : Stop position on the chromosome
        *value*       : Value for the feature
        *strand*      : Strand of the feature, '+' or '-'
        *phase*       : Phase of the feature if it's a CDS
        *attributes*  : A semicolon-delimited list of "field=data" strings.
                        For example,
                        'ID=FBtr0000123;Parent=FGgn0001;Name=transcript 1'
        
        *strvals*     :  By default, GFFFeature objects will have their
                         attributes typecast to integers, floats, etc.  However
                         if *strvals* is True, ALL attributes will be strings.
                         For example, with *strvals=False*, GFFFeature.start
                         and GFFFeature.end be integers (useful for downstream
                         work like ``featurelen = feature.stop -
                         feature.start``) but if *strvals=True* they will be
                         strings. 
        
                         Setting *strvals=True* will speed up parsing.

        *name*        : If provided, replaces the "chrom" kwarg.  Useful if
                        your GFF files are not "chromosome-centric".
        """
        if name is not None:
            if chrom is not None:
                raise ValueError, "specifying both chrom and name not supported"
            chrom = name
        if not strvals: # do typecasting
            self.chrom = chrom
            self.source = source
            self.featuretype = featuretype
            try:
                self.start = int(start)
            except (ValueError,TypeError):
                raise TypeError, 'start must be able to be converted to an integer'
            try:
                self.stop = int(stop)
            except (ValueError,TypeError):
                raise TypeError, 'stop must be able to be converted to an integer'
            if value is not None:
                try:
                    self.value = float(value)
                except (ValueError, TypeError):
                    self.value = None
            else:
                self.value = None
            self.strand = strand
            if phase is not None:
                try: 
                    self.phase = int(phase)
                except (ValueError, TypeError):
                    self.phase = None
            else:
                self.phase = None
          
        if strvals: # no typecasting, save everything as a string.
            self.chrom = chrom
            self.source = source
            self.featuretype = featuretype
            self.start = start
            self.stop = stop
            self.value = value
            self.strand = strand
            self.phase = phase

        self._parse_attributes(attributes)

    def _parse_attributes(self, attributes):  
        """
        Method to parse the attributes field of a line in a GFF file.
        """
        # keep track of these for later printing out.
        self._strattributes = attributes 

        if attributes is None:
            self._strattributes = ''
        # parse "attributes" field of the GFF line and insert them into an
        # Attributes object.
        self.attributes = GFFFeature.Attributes()
        if attributes is not None:
            items = attributes.split(';')
            for item in items:
                if len(item) == 0:
                    continue
                field,value = item.split('=')
                field = field.strip()
                values = value.split(',')
                values = [i.strip() for i in values]
                setattr(self.attributes, field,values)

                # Keep track inside the Attributes object of what you added to
                # it
                self.attributes._attrs.append(field)
        
    @property
    def id(self):
        try:
            return self.attributes.ID[0]
        except AttributeError:
            return None
    
    @property
    def chr(self):
        """Attribute *chr* now deprecated -- please use *chrom* instead"""
        return self.chrom

    @property
    def name(self):
        """Alias to chrom"""
        return self.chrom

    def add_attribute(self, attribute, value):
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
        if isinstance(value, basestring):
            value = [value]
        setattr(self.attributes, attribute, value)
        self.attributes._attrs.append(attribute)
        attr = ','.join(map(str,value))
        old_attrs = self._strattributes.split(';')
        old_attrs.append('%s=%s'%(attribute,attr))
        old_attrs = [i for i in old_attrs if len(i)>0]
        self._strattributes = ';'.join(old_attrs)+';'

        #self._strattributes += ';%s=%s' % (attribute, attr)

    def remove_attribute(self, attribute):
        """
        Delete attribute from this feature.  This method also removes the
        attribute from the equivalent GFF attributes string.  That is, if the
        starting attributes looked like this::

            ID=gene001;size=100kb;

        and you used remove_attribute('size'), then the new attributes string
        would look like::

            ID=gene001;
        """
        delattr(self.attributes, attribute)
        self.attributes._attrs.remove(attribute)
        ind1 = self._strattributes.find(attribute)
        ind2 = self._strattributes.find(';', ind1)
        self._strattributes = self._strattributes[:ind1] + self._strattributes[ind2:-1]

    def to_bed(self, fieldcount=3):
        """
        Returns the feature as a BED format line, with number of fields
        *fieldcount*.  Default is 3, which is chrom, start, stop.  Up to BED-6
        is supported.

        Note that a newline is added to the end of the string.  This allows for
        nice semantics like::
            
            >>> fout = open('genes.bed','w')
            >>> for gene in G.features_of_type('gene'):
            ...     fout.write(gene.to_bed(6))
            >>> fout.close()
        """
        attrs = ['chrom', 'start', 'stop', 'id', 'value', 'strand']
         
        fields = []
        if self.start is None:
            raise ValueError, 'feature start is None; need an integer to convert to BED'
        if self.stop is None:
            raise ValueError, 'feature start is None; need an integer to convert to BED'
        for i in range(fieldcount):
            if attrs[i] == 'start':
                fields.append(str(getattr(self, 'start')-1))
            else:
                fields.append(str(getattr(self, attrs[i])))
        return '\t'.join(fields)+'\n'

    def __repr__(self):
        return "%s %s '%s': %s:%s-%s (%s)" % (self.__class__.__name__, 
                                              self.featuretype,
                                              self.id,
                                              self.chrom,
                                              self.start,
                                              self.stop, 
                                              self.strand)

    def __len__(self):
        length = self.stop-self.start+1
        if length<1:
            raise ValueError, 'Zero- or negative length feature'
        return length
    
    def __eq__(self,other):
        """
        Test equality of two features based on their tostring() method, which 
        is the most uniform representation of a feature.
        """
        if not isinstance(other,self.__class__):
            raise ValueError, "cannot test equality to another object that's not the same class"
        if (self.chrom == other.chrom) &\
           (self.start == other.start) &\
           (self.stop == other.stop) &\
           (self.strand == other.strand) &\
           (self._strattributes == other._strattributes) &\
           (self.value == other.value) &\
           (self.id == other.id):
           return True
        else:
            return False

    def __ne__(self,other):
        return not self.__eq__(other)

    @property
    def TSS(self):
        """
        The transcription start site of the feature.  This is simply the
        strand-specific start of the feature.  Returns None if no strand
        specified.
        """
        if self.strand == '+':
            return self.start
        if self.strand == '-':
            return self.stop
        else:
            raise ValueError, 'TSS not defined for feature with strand=%s' % self.strand
    
    @property
    def midpoint(self):
        """
        Convenience accessor for getting the midpoint of a feature.
        """
        try:
            return self.start + (self.stop-self.start)/2
        except (TypeError,ValueError):
            raise ValueError, 'feature start and stop must be integers'
            

    def tostring(self):
        """
        Prints the GFF record suitable for writing to file (newline included).
        
        Since the string output is reconstructed based on the current contents
        of the GFFFeature, (attributes are reconstructed as well), this
        provides an easy means of editing GFF files, e.g.::
            
            >>> # shorten all CDS features by 10 bp, rename features to "*.short",
            >>> # and write only these shortened CDSs to file.
            >>> fout = open('out.gff')
            >>> for feature in GFFFile('in.gff'):
            ...     if feature.featuretype != 'CDS':
            ...         continue
            ...     feature.stop -= 10
            ...     feature.attributes.ID += '.short'
            ...     fout.write(feature.tostring())
            >>> fout.close()

        In the interest of speed, does not do any error-checking.
        """
        # Reconstruct the attributes field
        attributes = []
        for attr in self.attributes._attrs:
            values = getattr(self.attributes, attr)
            values = map(str, values)
            values = ','.join(values)
            attributes.append(attr + '='+values)
        attributes = ';'.join(attributes)
        if len(attributes) > 0:
            attributes += ';'

        items = [self.chrom, 
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
            print feature.chrom
            print feature.start
            print feature.stop
            print feature.featuretype
            print feature.desc
            print 'length:', feature.stop-feature.start
        """

    featureclass = GFFFeature

    def __init__(self, fname, strvals=False):
        """
        *fname*

            GFF filename.  Can be a string filename (.gz files detected via extension)
            or, if *fname* is not a string type, it can be a file-like object (sys.stdout, 
            already-open file, StringIO, etc).

        *strvals*

            By default, GFFFeature objects will have their attributes typecast
            to integers, floats, etc.  If *strvals* is True, ALL attributes
            will be strings.  For example, with *strvals=False*,
            GFFFeature.start and GFFFeature.end be integers (useful for
            downstream work) but if *strvals=True* they will be strings. 

            Setting *strvals=True* will speed up parsing.
        """
        if type(fname) is str:
            self.stringfn = True
            if os.path.splitext(fname)[-1] == '.gz':
                self.file = gzip.open(fname)
            else:
                self.file = open(fname)
        else:
            self.stringfn = False
            self.file = fname
        self.strvals = strvals

    def __iter__(self):
        """
        Yields a GFFFeature object for each line.
        """
        f = self.file
        for line in f:
            # You've reached the end of the GFF file; this represents the start
            # of the optional sequence section
            if line.startswith('>'):
                raise StopIteration
            line = line.rstrip()
            if line.startswith('#') or len(line) == 0:
                continue
            L = line.rstrip().split('\t')
            args = [None for i in range(10)]
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
    
    def __repr__(self):
        return 'gfffile object (file=%s)' % (self.file)

class GTFFeature(GFFFeature):
    """
    Class to represent a GTF feature and its annotations. Subclassed from GFFFeature.
    """
    def __init__(self, id=None, chrom=None, source=None, featuretype=None, start=None, stop=None,
                 value=None, strand=None, phase=None, attributes=None, name=None, strvals=False):
        
        GFFFeature.__init__(self,
                            chrom=chrom,
                            source=source,
                            featuretype=featuretype,
                            start=start,
                            stop=stop,
                            value=value,
                            strand=strand,
                            phase=phase,
                            attributes=attributes,
                            name=name,
                            strvals=strvals)
        if attributes is None:
            self._strattributes = ''
        else:
            self._strattributes = attributes
        if id is not None:
            self.add_attribute('ID', id)
    
    def tostring(self,gtfdb=None,include_id=True):
        """
        Prints the GTF record suitable for writing to file (newline
        included).
        
        In the interest of speed, does not do error-checking.

        AB000123    Twinscan     CDS    193817    194022    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2"; 
        """
        # 

        # Genes and transcripts are not explicitly written in GTF files.  So if gene or mRNA featuretype is asked for, you need to 
        # traverse the db and get the childrend for these objects.
        # 
        # Note that GTF files do not have multiple attribute values like GFF
        # files, so this means you'll have to output a line for each feature
        if self.featuretype == 'gene' or self.featuretype == 'mRNA':
            if gtfdb is None:
                raise ValueError, "Need to specify a GFFutils.GTFDB if you want to get all the GTF lines for this gene"

            if self.featuretype == 'gene':
                children = []
                for transcript in gtfdb.children(self,level=1):
                    for child in gtfdb.children(transcript,level=1):
                        child.remove_attribute('transcript_id')
                        child.add_attribute('transcript_id',transcript.id)
                        children.append(child)
            
            if self.featuretype == 'mRNA':
                children = []
                for child in gtfdb.children(self,level=1):
                    child.remove_attribute('transcript_id')
                    child.add_attribute('transcript_id',self.id)
                    children.append(child)


            lines = []
            children.sort(key=lambda x: (x.start,x.stop,x.featuretype))
            for child in children:
                if not include_id:
                    child.remove_attribute('ID')
                lines.append(child.tostring())
            return ''.join(lines)
                
         
        # Reconstruct the attributes field
        if not include_id:
            feature = copy.deepcopy(self)
        else:
            feature = self

        attributes = ''
        for attr in feature.attributes._attrs:
            values = getattr(feature.attributes, attr)
            if type(values) is list:
                values = ','.join(map(str, values))
            attributes += attr+' '+'"' + str(values)+'"; '

        items = [feature.chrom, 
                 feature.source,
                 feature.featuretype,
                 feature.start, 
                 feature.stop, 
                 feature.value, 
                 feature.strand, 
                 feature.phase,
                 attributes]
        
        printables = []
        for item in items:
            if item is None:
                printables.append('.')
            else:
                printables.append(str(item))
        return '\t'.join(printables).rstrip()+'\n'

    def add_attribute(self, attribute, value):
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
        setattr(self.attributes, attribute, value)
        self.attributes._attrs.append(attribute)
        msg = 'Multiple values for one attribute not currently supported for \
        GTF features'
        assert len(value) == 1, msg 
        
        old_attrs = self._strattributes.split(';')
        old_attrs.append('%s "%s"'%(attribute,value[0]))
        old_attrs = [i for i in old_attrs if len(i)>0]
        self._strattributes = ';'.join(old_attrs)+';'
  
        #self._strattributes += ';%s "%s"' % (attribute, value[0])

    def remove_attribute(self, attribute):
        """
        Delete attribute from this feature.
        """
        delattr(self.attributes, attribute)
        self.attributes._attrs.remove(attribute)
        ind1 = self._strattributes.find(attribute)
        ind2 = self._strattributes.find(';', ind1)
        self._strattributes = self._strattributes[:ind1] + self._strattributes[ind2:-1]

    def _parse_attributes(self, attributes):
        """
        Parse the attributes.  This is where GTF differs from GFF format.
        """
        # keep track of these for later printing out.
        self._strattributes = attributes 

        # parse "attributes" field of the GTF line and insert them into an
        # Attributes object.
        self.attributes = GTFFeature.Attributes()
        if attributes is not None:
            items = attributes.split(';')
            for item in items:
                if len(item) == 0:
                    continue
                field, value = item.strip().split()

                value = value.replace('"','')
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    pass
                setattr(self.attributes, field, value)

                # Keep track inside the Attributes object of what you added to
                # it
                self.attributes._attrs.append(field)
        
           

class GTFFile(GFFFile):
    """Iterator object, that moves through features in a GTF-format file.  A
    new GTFFeature object is created for each line.  Subclassed from GFFFile.
    
    Usage::

        for feature in GTFFile('a.bed'):
            print feature.chrom
            print feature.start
            print feature.stop
            print feature.featuretype
            print feature.desc
            print 'length:', feature.stop-feature.start
        """
    featureclass = GTFFeature

    def __repr__(self):
        return 'GTFFile object (file=%s)' % (self.file)

class Genome:
    """
    Wrapper class for quickly getting a sequence within a chromosome. Inspired by ERANGE.

    This class creates a memory-map of the fasta file.  It indexes where the newlines are
    in the file and records how many bits it is into the file.  Access is then
    pretty quick, since we can seek to that bit in the file rather quickly.
    
    For now, FASTA files must have their sequence all on one line. Use the
    fasta_seqs_to_oneline() function to make a file like this.

    Example usage::

        >>> g = Genome('dm3.fa')
        >>> nucleotides = g.sequence('chr2L',12000,13000)
    """
    def __init__(self, fn, debug=False):
        """
        *fn* is a FASTA-format file.
        """
        # Keep track of the starting position of each chromosome in self.startinds
        self.startinds = {}
        self.chromorder = []
        self.chromfiles = {}
        self.chromlens = {}
        self.namestarts = {}
        f = open(fn)
        m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        ind1 = 0
        m.seek(0)
        last_chrom = None
        while True:
            # Find the next occurrence of the fasta header delimiter.
            ind1 = m.find('>')
        
            # Assume that just before the ">" was a newline which indicated the
            # end of the last sequence. This will be -1 if we're on the first
            # chromosome, or -2 if there are no more ">" left in the file.
            end_of_last_seq = ind1-1
            

            # find() returns -1 if nothing found, so if no more "> then break
            # out of the loop.
            if ind1 == -1:
                # assume the end of the last sequence is the end of the file.
                # But subtract 1, because we assume the file ends with a
                # newline.
                end_of_last_seq = m.size()-1
                break

            # scoot up to the ">", and find the position of the newline at the end.
            m.seek(ind1)
            ind2 = m.find('\n')
            
            # The name of the chrom is from just after the ">" to just before
            # the newline (recall Python's half-open intervals, hence ind2
            # instead of ind2-1)
            chrom = m[ind1+1:ind2]

            # keep track of the order of chroms -- this will be used along with
            # the startinds to figure out the lengths of chromosomes.
            self.chromorder.append(chrom)
            
            if debug:
                print chrom, ind1, ind2
            
            # The sequence of this chromosome starts at ind2+1.
            self.startinds[chrom] = ind2+1

            # keep track of where the name starts
            self.namestarts[chrom] = ind1

            # Now that we're on the next chrom we can fill in the length of the
            # last one.
            if last_chrom is not None:
                self.chromlens[last_chrom] = end_of_last_seq - self.startinds[last_chrom]
            last_chrom = chrom

            # move to the newline.
            m.seek(ind2)

        
        # fill in the length of this last one.
        self.chromlens[last_chrom] = end_of_last_seq - self.startinds[last_chrom]

        
        ### TODO: the code below is a sketch of what to do about multimapping
        #   multiline FASTAs...but not working at the moment.
        
        if 0: 
            # count the newlines in each chrom; assume that line lengths are equal
            startinds = sorted(self.startinds.items(),key=lambda x: x[1])
            self.newlines = {}
            self.chromlens = {}
            for i in range(len(startinds)):
                chrom,startind = startinds[i]
                try:
                    nextchrom,nextstart = startinds[i+1]
                    nextnamestart = self.namestarts[nextchrom]
                    chromlen = nextnamestart-startind
                except IndexError:
                    nextnamestart = -1
                    chromlen = m.size()-startind 
                m.seek(startind)
                entire_chrom = m.read(nextnamestart-startind)
                self.newlines[chrom] = entire_chrom.count('\n')
                self.chromlens[chrom] = chromlen

        m.seek(0)
        self.mmap = m

    def sequence(self,chrom,start,stop,strand=None):
        """
        Returns the sequence for the position (*chrom*, *start*, *stop*)
        requested.

        If *strand* is '-', the reverse complement is returned.
        """
        if stop > self.chromlens[chrom]:
            raise ValueError,'stop position %s out of range for chrom %s of len %s)' % (stop,chrom,self.chromlens[chrom])

        # the chromosome start position
        i = self.startinds[chrom]
        
        # count the number of newlines between the chromosome start and the
        # position you want to get to.
        self.mmap.seek(0)

        start = i + start - 1
        self.mmap.seek(start)

        length = i + stop - start
        seq = self.mmap.read(length)

        if strand == '-':
            return seq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
        else:
            return seq
    
    def sequence_from_feature(self,feature):
        """
        Similar to self.sequence(), but accepts a GFFFeature or GTFFeature
        object for convenience.
        """
        return self.sequence(feature.chrom,feature.start,feature.stop,feature.strand)

    def splice_junctions(self,featureid,gffdb,padding=10):
        """
        Returns an iterator of spliced sequences for an mRNA, with *padding*
        additional bp on either side of the junction (which is itself 2bp).

        The *gffdb* needs to be specified so that the children exons can be
        queried.  If you want to get splices for an entire gene, just call this
        individually on each of that gene's mRNAs.
        """
        if type(featureid) is not str:
            featureid = featureid.id
        transcript = gffdb[featureid]
        exons = list(gffdb.children(featureid,featuretype='exon',level=1))

        exons.sort(key=lambda x: x.start)
        for i,exon in enumerate(exons):
            # we're on the last exon, which is not spliced to anything -- so
            # we're done.
            if i == len(exons)-1:
                break
            next_exon = exons[i+1]
            seq1 = self.sequence(chrom=exon.chrom,
                                 start=exon.stop-padding,
                                 stop=exon.stop,
                                 strand=exon.strand)
            seq2 = self.sequence(chrom=next_exon.chrom,
                                 start=next_exon.start,
                                 stop=next_exon.start+padding, 
                                 strand=next_exon.strand)
            yield seq1+seq2

            

    def spliced_transcript(self,featureid,gffdb,featuretype='exon'):
        """
        Returns the spliced sequence of an mRNA, *featureid*.  *featureid* can
        be either a string or a GFFFeature.

        The *gffdb* needs to be specified so that the children exons can be
        queried.
        
        *featuretype* is by default 'exon', which will give you the spliced
        transcript along with UTRs (assuming the UTRs are annotated as exons,
        which is the case in FlyBase GFF files).  If you'd prefer just the
        coding sequence, then use featuretype='CDS'.
        """
        if type(featureid) is not str:
            featureid = featureid.id

        # get the exons and make sure they're in order
        exons = list(gffdb.children(featureid,featuretype=featuretype,level=1))
        exons.sort(key=lambda x: x.start)

        seq = []
        for exon in exons:
            seq.append(self.sequence(chrom=exon.chrom,
                                     start=exon.start,
                                     stop=exon.stop,
                                     strand=exon.strand))
        return ''.join(seq)
            

class DBCreator(object):
    def __init__(self,dbfn):
        self.dbfn = dbfn
        self.Feature = self.__class__.featureclass
        conn = sqlite3.connect(dbfn)
        conn.text_factory = sqlite3.OptimizedUnicode
        self.conn = conn

    def drop_indexes(self):
        c = self.conn.cursor()
        c.execute('DROP INDEX IF EXISTS ids')
        c.execute('DROP INDEX IF EXISTS parentindex')
        c.execute('DROP INDEX IF EXISTS childindex')
        self.conn.commit()
    
    def init_tables(self):
        c = self.conn.cursor()
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
        self.conn.commit()

    def create_from_file(self,fn):
        if self.__class__.featureclass == GFFFeature:
            features = GFFFile(fn,strvals=False)
        if self.__class__.featureclass == GTFFeature:
            features = GTFFile(fn,strvals=True)
        self.init_tables()
        self.populate_from_features(features)
        self.update_relations()

class GFFDBCreator(DBCreator):
    featureclass = GFFFeature
    def __init__(self,dbfn):
        DBCreator.__init__(self,dbfn)

    def populate_from_features(self,features):
        c = self.conn.cursor()
        self.drop_indexes()
        for feature in features:
            if feature.id is None:
                new_id = '%s:%s:%s-%s' % (feature.featuretype, feature.chrom, feature.start, feature.stop)
                feature.add_attribute('ID',new_id)
            c.execute('''
                      INSERT OR IGNORE INTO features VALUES (?,?,?,?,?,?,?,?,?,?)
                      ''',(feature.id, 
                       feature.chrom,
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
        self.conn.commit()

    def update_relations(self):
        c = self.conn.cursor()
        c2 = self.conn.cursor()
        c3 = self.conn.cursor()
        c.execute('CREATE INDEX ids ON features (id)')
        c.execute('CREATE INDEX parentindex ON relations (parent)')
        c.execute('CREATE INDEX childindex ON relations (child)')
        self.conn.commit()

        c.execute('SELECT id FROM features')

        # create 2 more cursors so you can iterate over one while querying on
        # the other's iteration
        tmp = tempfile.mktemp()
        fout = open(tmp,'w')
        for parent in c:
            parent = parent[0] # first thing in the list is the ID
            
            # Here we get the first-level child from the initial import.   This
            # data was contained in the "Parent=" attribute of each GFF feature.
            c2.execute('SELECT child FROM relations WHERE parent = ? AND level=1', (parent,))
            for child in c2:
                child = child[0]
                c3.execute('SELECT child FROM relations WHERE parent = ? AND level=1', (child,))
                for grandchild in c3:
                    grandchild = grandchild[0]
                    fout.write('%s\t%s\n' % (parent,grandchild))
        fout.close()

        for line in open(tmp):
            parent,child = line.strip().split('\t')
            c.execute('INSERT OR IGNORE INTO relations VALUES (?,?,?)', (parent, child, 2))

        c.execute('drop index childindex')
        c.execute('drop index parentindex')
        c.execute('create index parentindex on relations (parent)')
        c.execute('create index childindex on relations (child)')
        c.execute('create index starts on features(start)')
        c.execute('create index stops on features(stop)')
        c.execute('create index startstrand on features(start, strand)')
        c.execute('create index stopstrand on features(stop,strand)')
        c.execute('create index featuretypes on features(featuretype)')

        self.conn.commit()
        os.unlink(tmp)

class GTFDBCreator(DBCreator):
    featureclass = GTFFeature
    def __init__(self,dbfn):
        DBCreator.__init__(self,dbfn)
    
    def populate_from_features(self,features):
        self.drop_indexes()
        c = self.conn.cursor()
        for feature in features:
            parent = feature.attributes.transcript_id
            grandparent = feature.attributes.gene_id
            
            # A database-specific ID to use
            ID = '%s:%s:%s-%s' % (feature.featuretype, feature.chrom, feature.start, feature.stop) 

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
                       feature.chrom,
                       feature.start, 
                       feature.stop,
                       feature.strand,
                       feature.featuretype,
                       feature.value,
                       feature.source,
                       feature.phase,
                       feature._strattributes))
        self.conn.commit()
    
    def update_relations(self):
        self.drop_indexes()
        c = self.conn.cursor()
        c2 = self.conn.cursor()
        c.execute('CREATE INDEX ids ON features (id)')
        c.execute('CREATE INDEX parentindex ON relations (parent)')
        c.execute('CREATE INDEX childindex ON relations (child)')

        tmp = tempfile.mktemp()
        fout = open(tmp,'w')
        c.execute("SELECT DISTINCT parent FROM relations")
        for parent in c:
            parent = parent[0]
            c2.execute("""
                       SELECT min(start), max(stop), level, strand, chrom FROM features 
                       JOIN relations ON
                       features.ID = relations.child
                       WHERE
                       parent = ? 
                       AND featuretype == "exon"
                       """, (parent,))
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

        self.conn.commit()
        os.remove(tmp)
        
        c.execute('DROP INDEX ids')
        c.execute('CREATE INDEX ids ON features (id)')
        self.conn.commit()


def fasta_seqs_to_oneline(infn, outfn):
    """
    Converts a typical FASTA file *infn*, with many lines per sequence, into
    one with a single long line for each sequence and saves this as *outfn*.
    This new file can then be used with the Genome class for fast indexing into
    the genome.
    """
    fout = open(outfn,'w')
    nchrom = 0
    for line in open(infn):
        if line.startswith('>'):
            if nchrom == 0:
                fout.write(line)
            else:
                fout.write('\n')
                fout.write(line)
            nchrom += 1
            continue
        else:
            fout.write(line.rstrip())
    # write a newline at the very end.
    fout.write('\n')
    fout.close()

def inspect_featuretypes(gfffn):
    """
    Returns a list of unique featuretypes in the GFF file, *gfffn*.  Useful for
    when you're about to clean up a GFF file and you want to know which
    featuretypes to exclude in the cleaned GFF.
    """
    featuretypes = set()
    for line in open(gfffn):
        L = line.split('\t')
        try:
            featuretypes = featuretypes.union([L[2]])
        except IndexError:
            continue
    return list(featuretypes)
        

def clean_gff(gfffn,newfn=None,addchr=False,featuretypes_to_remove=None,sanity_check=True):
    """
    Helps prepare a GFF file *gfffn* for import into a database. The new,
    cleaned file will be saved as *newfn* (by default, *newfn* will be *gfffn*
    plus a ".cleaned" extension).
    
    Use the inspect_featuretypes() function in this module to determine what
    featuretypes are contained in the GFF; then you can filter out those that
    are not interesting by including them in the list of
    *featuretypes_to_remove*

    If *addchr* is True, the string "chr" will be prepended to the chromosome
    name.

    If *sanity_check* is True, only features that pass a simple check will be
    output.  Currently the only sanity check is that coordinates are not
    negative and that feature starts are not greater than feature stop coords.

    Also, some GFF files have FASTA sequence at the end.  This function will
    remove the sequence as well.
    """
    if newfn is None:
        newfn = gfffn+'.cleaned'
    if featuretypes_to_remove is None:
        featuretypes_to_remove = []
    
    fout = open(newfn,'w')

    for line in open(gfffn):
        if line.startswith('>'):
            break
        L = line.split('\t')
    
        try:
            if L[2] in featuretypes_to_remove:
                continue
        except IndexError:
            continue
        
        if sanity_check:
            start = int(L[3])
            stop = int(L[4])
            if start<0 or stop<0 or (start>stop):
                continue
        
        if addchr:
            fout.write('chr'+line)
        else:
            fout.write(line)
    
    fout.close()

def create_gffdb(gfffn, dbfn, verbose=True):
    """
    Reads in a GFF3 file (with autodetection and support for *.gz files, based
    on extension) and constructs a database for use with downstream analysis.
    See the GFFDB class in particular for an interface to this database. 

    This takes 2-3 min on a 100 MB GFF file for D. melanogaster.  This is a
    one-shot time investment, since once it's created using the database is
    quite fast.

    You may want to try the clean_gff() function on your GFF file, which can
    add a 'chr' to the beginning of chromosome names (useful for FlyBase GFFs)
    or remove featuretypes from the GFF that you're not interested in.  This
    latter part can considerably speed up the database creation and subsequent
    use.

    The database requires features to have unique names (the ID attribute).  If
    there is no "ID" attribute, then the feature will be named according to its
    featuretype and an incrementing counter and with an "unnamed" prefix.  For
    example, the first exon in the file to not have an ID attribute will be
    called "unnamed_exon:1".

    Note that these labels are not stable across multiple GFF files!  Adding
    your own feature IDs directly to the GFF file is the most robust solution.


    IMPLEMENTATION:

    There are two tables in the database. 
    
    The "features" table contains the fields of the GFF line -- chrom, start,
    stop, strand, featuretype, value, source, phase, attributes.  It also
    contains the "id" field, which is a unique ID for each feature.
    
    Schema of the database:

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
        CREATE INDEX childindex on relations (child);
        CREATE INDEX featuretypes on features(featuretype);
        CREATE INDEX ids ON features (id);
        CREATE INDEX parentindex on relations (parent);
        CREATE INDEX starts on features(start);
        CREATE INDEX startstrand on features(start, strand);
        CREATE INDEX stops on features(stop);
        CREATE INDEX stopstrand on features(stop,strand);
    """

    # Calculate lines so you can display the percent complete
    if os.path.splitext(gfffn)[-1] == '.gz':
        f = gzip.open(gfffn)
    else:
        f = open(gfffn)
    nlines = 0.0
    for line in f:
        if line.startswith('>'):
            break
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
            if verbose:
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
                  INSERT OR IGNORE INTO features VALUES (?,?,?,?,?,?,?,?,?,?)
                  ''',(feature.id, 
                   feature.chrom,
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
    if verbose:
        print '(%ds)' % (t1-t0)
    t0 = time.time()

    # Index the database to speed up the second-level queries.
    if verbose:
        print 'creating indexes: 0%',
        sys.stdout.flush()
    c.execute('CREATE INDEX ids ON features (id)')
    if verbose:
        print '\rcreating indexes: 33%',
        sys.stdout.flush()
    c.execute('create index parentindex on relations (parent)')
    if verbose:
        print '\rcreating indexes: 66%',
        sys.stdout.flush()
    c.execute('create index childindex on relations (child)')
    if verbose:
        print '\rcreating indexes: 100%',
        sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    if verbose:
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
            if verbose:
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
    if verbose:
        print '(%ds)' % (t1-t0)

    t0 = time.time()
    counter = 0
    last_perc = 0
    grandchild_count = float(grandchild_count)
    for line in open(tmp):
        counter += 1
        perc_done = int(counter / grandchild_count * 100)
        if perc_done != last_perc:
            if verbose:
                print '\rimporting grandchildren: %s%%' % perc_done,
                sys.stdout.flush()
        last_perc = perc_done
        parent,child = line.strip().split('\t')
        c.execute('insert or ignore into relations values (?,?,?)', (parent, child, 2))
    t1 = time.time()
    if verbose:
        print '(%ds)' % (t1-t0)
    t0 = time.time()
    if verbose:
        print 're-creating indexes',
    c.execute('drop index childindex')
    c.execute('drop index parentindex')
    c.execute('create index parentindex on relations (parent)')
    if verbose:
        print '\rre-creating indexes: 50%',
        sys.stdout.flush()
    c.execute('create index childindex on relations (child)')
    c.execute('create index starts on features(start)')
    c.execute('create index stops on features(stop)')
    c.execute('create index startstrand on features(start, strand)')
    c.execute('create index stopstrand on features(stop,strand)')
    c.execute('create index featuretypes on features(featuretype)')
    if verbose:
        print '\rre-creating indexes: 100%',
        sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    if verbose:
        print '(%ds)' % (t1-t0)
    os.remove(tmp)
        

def create_gtfdb(gtffn, dbfn, verbose=True):
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
            if verbose:
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
        ID = '%s:%s:%s-%s' % (feature.featuretype, feature.chrom, feature.start, feature.stop)
        

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
                   feature.chrom,
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
    if verbose:
        print '(%ds)' % (t1-t0)

    # Now that the entries in the GTF file have been added, we can index
    # the database to speed up the second-level queries we need to do to
    # identify genes.
    t0 = time.time()
    if verbose:
        print 'creating indexes: 0%',
        sys.stdout.flush()
    c.execute('DROP INDEX IF EXISTS ids')
    c.execute('CREATE INDEX ids ON features (id)')
    if verbose:
        print '\rcreating indexes: 33%',
        sys.stdout.flush()
    c.execute('DROP INDEX IF EXISTS parentindex')
    c.execute('create index parentindex on relations (parent)')
    if verbose:
        print '\rcreating indexes: 66%',
        sys.stdout.flush()
    c.execute('drop index if exists childindex')
    c.execute('CREATE INDEX childindex ON relations (child)')
    if verbose:
        print '\rcreating indexes: 100%',
        sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    if verbose:
        print '(%ds)' % (t1-t0)

    
    # GTF files only describe exons, so genes and transcripts are described
    # implicitly.  With all the exons imported into the database, we can
    # now do queries to find the limits of genes and transcripts . . . 

    # Get an idea of how many parents we have to get through
    t0 = time.time()
    c.execute("SELECT COUNT(DISTINCT parent) FROM relations")
    nparents = float(c.fetchone()[0])
    
    # These results will include only transcripts and genes, since exons and
    # CDSs aren't parents of anything and so are not in this table as parents
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
            if verbose:
                print '\rquerying database for start/stop positions of genes and transcripts: %s%%'%perc_done,
                sys.stdout.flush()
        last_perc = perc_done

        # Find the start and stop limits of this parent's children
        parent = parent[0]

        # parent may be None if there's something funky in the GTF file.
        if parent is None:
            continue


        c2.execute("""
                   select min(start), max(stop), level, strand, chrom FROM features 
                   JOIN relations ON
                   features.ID = relations.child
                   WHERE
                   parent = ? 
                   AND featuretype == "exon"
                   """, (parent,))
        
        # For testing to make sure you only get one level back.
        #child_limits = c2.fetchall()
        #assert len(child_limits) == 1
        #child_limits = child_limits[0]

        child_limits = c2.fetchone()
        start,end,level,strand, chrom = child_limits
        
        #print parent, start, end, level

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
    if verbose:
        print '(%ds)' % (t1-t0)

    # Re-create indexes since we just added all those parents. Prob only
    # need to re-create the first one though.
    t0 = time.time()
    if verbose:
        print 're-creating indexes: 0%',
        sys.stdout.flush()
    c.execute('DROP INDEX ids')
    c.execute('CREATE INDEX ids ON features (id)')
    if verbose:
        print '\rre-creating indexes: 33%',
        sys.stdout.flush()
    conn.commit()
    t1 = time.time()
    if verbose:
        print '(%ds)' % (t1-t0)
    conn.commit()
    del feature_counts

class GFFDB:
    featureclass = GFFFeature
    add_id = ''
    def __init__(self,db_fn):
        """
        Class to interface with a database created by the create_gffdb()
        function.
        
        *db_fn* is the filename of the database to use.  See the help for
        create_gffdb() for much more information on the creation and design
        ideas behind the database.
        """
        self.db_fn = db_fn
        self.conn = sqlite3.connect(db_fn)
        self.conn.text_factory = str

    def __getitem__(self,id):
        if isinstance(id, self.__class__.featureclass):
            id = id.id
        c = self.conn.cursor()
        c.execute('''
                  SELECT %s chrom, source, featuretype, start, stop, value,
                  strand, phase, attributes from features where id = ?
                  ''' % self.__class__.add_id, (id,))
        results = c.fetchall()
        if not len(results) == 1:
            raise FeatureNotFoundError(id)
        return self.__class__.featureclass(*results[0])

    def attribute_search(self,text,featuretype='gene'):
        """
        Looks for *text* within the "attributes" field of features of
        *featuretype* ('gene' by default).  Useful for looking up genes based
        on symbol or name rather than accession number.  Returns a list of
        matching GFFFeatures.

        Uses SQL's LIKE operator, which is case-insensitive.
        """
        text = '%'+text+'%'
        c = self.conn.cursor()
        c.execute('''
                  SELECT %s chrom, source, featuretype, start, stop, value,
                  strand, phase, attributes from features where attributes
                  like "%s" and featuretype = ?
                  ''' % (self.__class__.add_id, text),(featuretype,))
        results = []
        for result in c.fetchall():
            results.append(self.__class__.featureclass(*result))
        return results

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
        SQL queries . . . or just call this method iteratively.
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

    def exonic_bp(self,id,ignore_strand=False):
        """
        Merges all exons of the gene *id* and sums the total bp.
        *ignore_strand* is passed to merge_features(). Useful for calculating
        RPKM for full gene models.
        """
        if isinstance(id, self.__class__.featureclass):
            id = id.id
        exons = self.children(id, featuretype='exon', level=2)
        exons = list(exons)
        if len(exons) == 0:
            return 0
        merged_exons = self.merge_features(exons,ignore_strand)
        total_exon_bp = 0
        for exon in merged_exons:
            total_exon_bp += len(exon)
        return total_exon_bp

    def closest_feature(self,chrom,pos,strand=None,featuretype=None, ignore=None):
        """
        Returns the closest feature and the distance to the start position of
        that feature from *pos*.  It's up the caller to determine things like
        upstream/downstream, TSS, TTS, midpoints, etc.

        *ignore* is either None (default) or a list of IDs that you would like
        to ignore.
        
        If you are providing the *chrom*, *pos* of a gene's TSS, then you will
        have to add that gene's ID to the *ignore* list so you don't get that
        same gene returned as the closest (which it is).  
        
        In this case, you'll probably want to provide *featuretype* too (e.g.,
        featuretype='gene') -- otherwise you'll also have to add the gene's
        mRNAs, exons, introns, CDSs, and any other related features to the
        *ignore* list.

        The *ignore* list also makes it possible to get the second-closest
        feature to a coord by calling this method to get the closest feature,
        then adding that closest feature's ID to the *ignore* list and calling
        this method a second time.

        Note that providing a *featuretype* and *strand* will speed things up
        considerably.

        To illustrate, here is a somewhat contrived example usage to get
        closest upstream TSS on the plus strand to the genomic position
        "chr2L:5000000" (and stop looking if there's nothing within 1Mb)::

            >>> chrom = 'chr2L'
            >>> pos = 5000000
            >>> closest_pos = pos+1
            >>> closest_feature = None
            >>> strand = None
            >>> ignore = []
            >>> while (closest_pos > pos) or (strand !='-'):
            ...     dist,closest_feature = G.closest_feature(chrom=chrom,pos=pos,ignore=ignore,featuretype='gene')
            ...     strand = closest_feature.strand
            ...     ignore.append(closest_feature.id)
            ...     closest_pos = closest_feature.TSS
            ...     if dist > 1e6:
            ...         closest_feature = None
            ...         break


        """
        # e.g., AND id != FBgn0001 AND id != FBgn0002
        ignore_clause = ''
        if ignore is not None:
            if type(ignore) is str:
                ignore = [ignore]
            for i in ignore:
                ignore_clause += ' AND id != "%s" ' % i
        
        strand_clause = ''
        if strand is not None:
            strand_clause = ' AND strand="%s" ' % strand
        
        featuretype_clause = ''
        if featuretype is not None:
            featuretype_clause = ' AND featuretype="%s" ' % featuretype
        
        c = self.conn.cursor()
        
        c.execute('''
        SELECT abs(%(pos)s-start) as AAA,id FROM features 
        WHERE
        chrom = ?
        %(ignore_clause)s
        %(strand_clause)s
        %(featuretype_clause)s
        ORDER BY AAA LIMIT 1
        ''' % locals(),(chrom,))
        closest_start = c.fetchone()
        
        c.execute('''
        SELECT abs(%(pos)s-stop) as AAA,id FROM features 
        WHERE
        chrom = ?
        %(ignore_clause)s
        %(strand_clause)s
        %(featuretype_clause)s
        ORDER BY AAA LIMIT 1
        ''' % locals(),(chrom,))
        closest_stop = c.fetchone()

        candidates = [closest_start,closest_stop]
        candidates = [i for i in candidates if i is not None] 
        if len(candidates) == 0:
            return None,None

        candidates.sort(key=lambda x: x[0])
        dist,feature_id = candidates[0]
        return dist, self[feature_id]

    def closest_TSS(self, chrom, pos, featuretype='gene', strand=None, ignore=None, direction=None):
        raise DeprecationWarning, "closest_TSS() is now deprecated -- please use closest_feature() instead, which relies on the caller for more complex manipulation"

    def overlapping_features(self, chrom, start, stop, featuretype=None, strand=None, completely_within=False):
        """
        Returns an iterator of features of type *featuretype* that overlap the
        coords provided. If *featuretype* is None (default), then all features
        will be returned.

        If *strand* is not None, then only features of the strand specifed will be returned.

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
        Returns an iterator of merged features (overlapping features are
        combined into one, spanning from the lowest start to the highest stop
        position)

        *features* is an iterable of features that you want to merge together.

        *features* will be converted into a list because this method needs to
        sort the features.

        If *ignore_strand* is True, strand will be forced to all '+' and no checking
        will be done on the input

        If *ignore_strand* is False (default), then *features* must be all on
        the same strand.  
        
        The new *featuretype* will be a joining of the various featuretypes
        that were provided.  If all the features were 'exon', then the new,
        merged feature will have a featuretype of 'merged_exon'.  If there were
        some introns in *features*, then the new merged feature will have a
        featuretype of 'merged_exon_intron'.
        
        Note that merged features are not saved to the database -- they are
        created on-the-fly and only exist in memory.
        """
        # If it quacks like an iterator...then turn it into a list.
        if hasattr(features,'next'):
            features = list(features)
        
        # Quick function to sort by start position
        def start_pos(x):
            return x.start
        
        # Sort 'em by start position
        features.sort(key=start_pos)
        
        # Either set all strands to '+' or check for strand-consistency.
        if ignore_strand:
            strand = '+'
        else:
            strands = [i.strand for i in features]
            if len(set(strands))> 1:
                raise ValueError, 'Specify ignore_strand=True to force merging of multiple strands'
            strand = strands[0]
        
        # Sanity check to make sure all features are from the same chromosome.
        chroms = [i.chrom for i in features]
        if len(set(chroms)) > 1:
            raise NotImplementedError, 'Merging multiple chromosomes not implemented'
        chrom = chroms[0]

        # If they were all exons, merged objects will have featuretype
        # 'merged_exon'.  If features included both exon and CDS, will have
        # featuretype 'merged_exon_CDS'.
        featuretypes = list(set([i.featuretype for i in features]))
        featuretypes.sort()
        featuretypes = map(str,featuretypes)
        featuretypes = '_'.join(featuretypes)
        featuretype = 'merged_%s' % featuretypes

        # To start, we create a merged feature of just the first feature.
        current_merged_start = features[0].start
        current_merged_stop = features[0].stop

        # We don't need to check the first one, so start at feature #2.
        for feature in features[1:]:
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
                # The start position is outside the merged feature, so we're done with 
                # the current merged feature.  Prepare for output...
                Feature = feature.__class__
                merged_feature = Feature(chrom=feature.chrom,
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
        merged_feature = Feature(chrom=feature.chrom,
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

        So if *features* contains the exons of a transcript, this method will
        return the introns.
        
        This is a purposefully naive method that does NOT do sorting or merging
        -- it simply returns the intervening space between the features provided.

        If the features returned are not sorted, you may get overlapping
        results.

        Example for getting all the 'non-exonic' space on the plus strand of
        chrX::
            
            >>> # merge_features needs a single chrom and single strand
            >>> exons = G.features_of_type('exon',chrom='chrX',strand='+')
            >>> merged_exons = G.merge_features(exons)
            >>> non_exonic = G.interfeatures(merged_exons)

        Example for getting the introns of a transcript::

            >>> transcript = G.random_feature('mRNA')
            >>> exons = G.children(transcript,'exon')
            >>> introns = G.interfeatures(exons)
        """
        for i,feature in enumerate(features):
            if i == 0:
                interfeature_start = feature.stop
                last_feature = feature
                continue
            interfeature_stop = feature.start
            featuretype = 'inter_%s_%s' % (last_feature.featuretype, feature.featuretype)
            assert last_feature.strand == feature.strand
            assert last_feature.chrom == feature.chrom
            strand = last_feature.strand
            chrom = last_feature.chrom

            # Shrink
            interfeature_start += 1
            interfeature_stop -= 1

            yield self.featureclass(chrom=chrom,
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
 
    def children(self,id,level=1,featuretype=None):
        """
        Returns an iterator of the children *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified.
        """

        if level > 2:
            raise NotImplementedError, 'Levels > 2 not supported yet.'

        if isinstance(id, self.__class__.featureclass):
            id = id.id

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
        """
        Returns an iterator of the parents *level* levels away.  Immediate
        children are level=1; "grandchildren" are level=2.  The child
        transcripts of a gene will be at level=1; the child exons of the same
        gene are at level=2.  Returns an error if there are not enough levels
        of children for the level specified."""

        if level > 2:
            raise NotImplementedError, 'Levels > 2 not supported yet.'

        if isinstance(id, self.__class__.featureclass):
            id = id.id

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

            cdsStart and cdsEnd are the min and max positions, respectively, of all CDSs in the _gene_.
            So a particular isoform's CDS may not actually start on the gene-wide minimum CDS position.

            Assumes that there was an attribute in the GFF file called 'Name'.  This will then be
            found in the gene.attributes.Name attribute.
        """
        
        geneID = gene.id
        for transcript in self.children(geneID,level=1,featuretype='mRNA'):
            
            txStart = transcript.start
            txEnd = transcript.stop
            exons = []
            cdss = []
            exon_count = 0
            for i in self.children(transcript,1):
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

            exonStarts = ','.join([str(i.start) for i in exons])+','
            exonEnds = ','.join([str(i.stop) for i in exons])+','
            try:
                name = transcript.attributes.Name[0]
            except AttributeError:
                name = transcript.id
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (name, 
                                                                     gene.id,
                                                                     gene.chrom, 
                                                                     gene.strand, 
                                                                     gene.start, 
                                                                     gene.stop,
                                                                     cdsStart,
                                                                     cdsEnd,
                                                                     len(exons),
                                                                     exonStarts,
                                                                     exonEnds)
            yield line

    def count_features_of_type(self,featuretype):
        """
        More memory efficient than the functionally-equivalent::

            >>> len(list(self.features_of_type(featuretype)))
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
        include.  The returned feature will have a legnth of *dist* + 1, since
        the TSS will be included as well.
        
        TSS is considered the feature start position if the feature is on the
        plus strand, or the feature stop position if on the minus strand. The
        application of *dist* to getting a promoter can be changed by the other
        kwargs below:

        If *bidirectional* is True (default), then include *dist* bp on BOTH
        sides of the TSS.  Otherwise, only include the region *dist* bp 5' of
        the TSS.
        
        If *truncate_at_next_feature* (None by default) is a featuretype (e.g.,
        'gene' or 'mRNA') then first try and go up to *dist* bp upstream from
        feature start.  If there is another feature of the specified
        featuretype within *dist* upstream on either strand, then the promoter
        region returned will be truncated to stop at this feature.  This is
        useful if you don't want your promoter regions to overlap with another
        gene.
        
        If *truncate_at_next_feature* is not None AND *bidirectional* is True,
        then the downstream distance will not be truncated.

        Note that without supplying the chromosome limits, it's impossible to
        know whether the returned promoter runs off the end of the chromosome.
        The lower boundary is accounted for (i.e., no negative coordinates) but
        checking the upper boundary is the responsibility of the calling
        function.
        """
        if not isinstance(id, self.__class__.featureclass):
            feature = self[id]
        else:
            feature = id

        c = self.conn.cursor()
        
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
                                       feature.chrom, 
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
                                      feature.chrom, 
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

        # Create a new GFFFeature object (or GTFFeature) to return
        if self.__class__.featureclass == GTFFeature:
            promoter = self.__class__.featureclass(id='promoter:'+id,
                                          chrom=feature.chrom,
                                          source='imputed',
                                          featuretype='promoter',
                                          start=start,
                                          stop=stop,
                                          strand=feature.strand,
                                          value=None,
                                          phase=None,
                                          attributes=None)
     
        else:
            promoter = self.__class__.featureclass(chrom=feature.chrom,
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
        featuretype_subclause = ''
        if featuretype is not None:
            featuretype_clause = 'featuretype = "%s" AND ' % featuretype
            featuretype_subclause = 'WHERE featuretype = "%s"' % featuretype
            featuretype_subclause = ''
        c.execute('''
        SELECT %s chrom, source, featuretype, start, stop, value, strand, phase, attributes 
        FROM features
        WHERE 
        %s
        rowid >= abs(random()) %% (SELECT count() FROM features %s) LIMIT 1
        ''' % (self.__class__.add_id, featuretype_clause, featuretype_subclause))
        results = c.fetchone()
        return self.__class__.featureclass(*results)
        

    def coding_genes(self):
        """
        Returns an iterator of genes that contain CDSs. Useful for if you want
        to exclude tRNA, various ncRNAs, etc, since they are also annotated
        with featuretype "gene" and contain 'exon' features (at least in
        FlyBase GFFs)
        """
        for g in self.features_of_type('gene'):
            coding = False
            for grandchild in self.children(g.id, level=2):
                if grandchild.featuretype == 'CDS':
                    coding = True
            if coding:
                yield g
    
    def n_gene_isoforms(self, geneID):
        """
        Returns the number of isoforms that this gene has.
        """
        if type(geneID) is not str:
            geneID = geneID.id

        n = 0
        for i in self.children(geneID, level=1, featuretype='mRNA'):
            n += 1
        return n

    def n_exon_isoforms(self, exonID):
        """
        Returns the number of isoforms that this exon is found in.
        """

        if type(exonID) is not str:
            exonID = exonID.id

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

        if type(exonID) is not str:
            exonID = exonID.id

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
    
    def __init__(self, db_fn):
        GFFDB.__init__(self,db_fn)

    def UTRs(self, transcript_id):
        """
        Returns 5' and 3' UTRs for a transcript
        """
        raise NotImplementedError, 'Still working on the tests for this method.'
        transcript = self[transcript_id]
        
        exons = list(self.children(transcript_id, level=1, featuretype='exon'))
        exons.sort(key=lambda x: x.start)
        cdss = list(self.children(transcript_id, level=1, featuretype='CDS'))
        cdss.sort(key=lambda x: x.start)

        UTRs = []
        if transcript.strand == '+':
            first_cds = cdss[0]
            last_cds = cdss[-1]
            for exon in exons:
                # this can happen with spliced UTRs:
                if exon.stop < first_cds.start:
                    UTR = self.__class__.featureclass(chrom=exon.chrom,
                                                      start=exon.start,
                                                      stop=exon.stop,
                                                      strand=exon.strand,
                                                      featuretype='five_prime_UTR')
                    UTRs.append(UTR)
                
                # here's the canonical case, where a CDS is a subset of an exon.
                if exon.start < first_cds.start < exon.stop:
                    UTR = self.__class__.featureclass(chrom=exon.chrom,
                                                      start=exon.start,
                                                      stop=first_cds.start-1,
                                                      strand=exon.strand,
                                                      featuretype='five_prime_UTR')
                    UTRs.append(UTR)

                if exon.start < last_cds.stop < exon.stop:
                    UTR = self.__class__.featureclass(chrom=exon.chrom,
                                                      start=last_cds.stop+1,
                                                      stop=exon.stop,
                                                      strand=exon.strand,
                                                      featuretype='three_prime_UTR')
                    UTRs.append(UTR)
        
                if exon.stop > last_cds.stop:
                    UTR = self.__class__.featureclass(chrom=exon.chrom,
                                                      start=exon.start,
                                                      stop=exon.stop,
                                                      strand=exon.strand,
                                                      featuretype='three_prime_UTR')

        return UTRs
