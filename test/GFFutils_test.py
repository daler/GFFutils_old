import GFFutils
import os
import sqlite3
import nose.tools as nt
import difflib
import pprint
import copy


testdbfn = 'testing.db'

def EXPECTED_DATA():
    # list the children and their expected first-order parents for the GFF test file.
    GFF_parent_check_level_1 = {'FBtr0300690':['FBgn0031208'],
                                'FBtr0300689':['FBgn0031208'],
                                'exon:chr2L:7529-8116':['FBtr0300689','FBtr0300690'],
                                'five_prime_UTR_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                                'CDS_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                                'intron_FBgn0031208:1_FBgn0031208:2':['FBtr0300690'],
                                'intron_FBgn0031208:1_FBgn0031208:3':['FBtr0300689'],
                                'FBgn0031208:3':['FBtr0300689'],
                                'CDS_FBgn0031208:3_737':['FBtr0300689'],
                                'CDS_FBgn0031208:2_737':['FBtr0300690'],
                                'FBgn0031208:2':['FBtr0300690'],
                                'intron_FBgn0031208:2_FBgn0031208:4':['FBtr0300690'],
                                'three_prime_UTR_FBgn0031208:3_737':['FBtr0300689'],
                                'FBgn0031208:4':['FBtr0300690'],
                                'CDS_FBgn0031208:4_737':['FBtr0300690'],
                                'three_prime_UTR_FBgn0031208:4_737':['FBtr0300690'],
                               }

    # and second-level . . . they should all be grandparents of the same gene.
    GFF_parent_check_level_2 = {
                                'exon:chr2L:7529-8116':['FBgn0031208'],
                                'five_prime_UTR_FBgn0031208:1_737':['FBgn0031208'],
                                'CDS_FBgn0031208:1_737':['FBgn0031208'],
                                'intron_FBgn0031208:1_FBgn0031208:2':['FBgn0031208'],
                                'intron_FBgn0031208:1_FBgn0031208:3':['FBgn0031208'],
                                'FBgn0031208:3':['FBgn0031208'],
                                'CDS_FBgn0031208:3_737':['FBgn0031208'],
                                'CDS_FBgn0031208:2_737':['FBgn0031208'],
                                'FBgn0031208:2':['FBgn0031208'],
                                'intron_FBgn0031208:2_FBgn0031208:4':['FBgn0031208'],
                                'three_prime_UTR_FBgn0031208:3_737':['FBgn0031208'],
                                'FBgn0031208:4':['FBgn0031208'],
                                'CDS_FBgn0031208:4_737':['FBgn0031208'],
                                'three_prime_UTR_FBgn0031208:4_737':['FBgn0031208'],
                               }

    # Same thing for GTF test file . . .
    GTF_parent_check_level_1 = {                                
                                'exon:chr2L:7529-8116':['FBtr0300689','FBtr0300690'],
                                'exon:chr2L:8193-9484':['FBtr0300689'],
                                'exon:chr2L:8193-8589':['FBtr0300690'],
                                'exon:chr2L:8668-9484':['FBtr0300690'],
                                'exon:chr2L:10000-11000':['transcript_Fk_gene_1'],
                                'exon:chr2L:11500-12500':['transcript_Fk_gene_2'],
                                'CDS:chr2L:7680-8116':['FBtr0300689','FBtr0300690'],
                                'CDS:chr2L:8193-8610':['FBtr0300689'],
                                'CDS:chr2L:8193-8589':['FBtr0300690'],
                                'CDS:chr2L:8668-9276':['FBtr0300690'],
                                'CDS:chr2L:10000-11000':['transcript_Fk_gene_1'],
                                'FBtr0300689':['FBgn0031208'],
                                'FBtr0300690':['FBgn0031208'],
                                'transcript_Fk_gene_1':['Fk_gene_1'],
                                'transcript_Fk_gene_2':['Fk_gene_2'],
                                'start_codon:chr2L:7680-7682':['FBtr0300689','FBtr0300690'],
                                'start_codon:chr2L:10000-11002':['transcript_Fk_gene_1'],
                                'stop_codon:chr2L:8611-8613':['FBtr0300689'],
                                'stop_codon:chr2L:9277-9279':['FBtr0300690'],
                                'stop_codon:chr2L:11001-11003':['transcript_Fk_gene_1'],
                                }


    GTF_parent_check_level_2 = {
                                'exon:chr2L:7529-8116':['FBgn0031208'],
                                'exon:chr2L:8193-9484':['FBgn0031208'],
                                'exon:chr2L:8193-8589':['FBgn0031208'],
                                'exon:chr2L:8668-9484':['FBgn0031208'],
                                'exon:chr2L:10000-11000':['Fk_gene_1'],
                                'exon:chr2L:11500-12500':['Fk_gene_2'],
                                'CDS:chr2L:7680-8116':['FBgn0031208'],
                                'CDS:chr2L:8193-8610':['FBgn0031208'],
                                'CDS:chr2L:8193-8589':['FBgn0031208'],
                                'CDS:chr2L:8668-9276':['FBgn0031208'],
                                'CDS:chr2L:10000-11000':['Fk_gene_1'],
                                'FBtr0300689':['FBgn0031208'],
                                'FBtr0300690':['FBgn0031208'],
                                'transcript_Fk_gene_1':['Fk_gene_1'],
                                'transcript_Fk_gene_2':['Fk_gene_2'],
                                'start_codon:chr2L:7680-7682':['FBgn0031208'],
                                'start_codon:chr2L:10000-11002':['Fk_gene_1'],
                                'stop_codon:chr2L:8611-8613':['FBgn0031208'],
                                'stop_codon:chr2L:9277-9279':['FBgn0031208'],
                                'stop_codon:chr2L:11001-11003':['Fk_gene_1'],
                               }

    expected_feature_counts = {
                'GFF':{'gene':3,
                       'mRNA':4,
                       'exon':6,
                       'CDS':5,
                       'five_prime_UTR':1,
                       'intron':3,
                       'pcr_product':1,
                       'protein':2,
                       'three_prime_UTR':2},
                'GTF':{'gene':3,
                       'mRNA':4,
                       'CDS':5,
                       'exon':6,
                       'start_codon':2,
                       'stop_codon':3}
                }

    expected_features = {'GFF':['gene',
                                'mRNA',
                                'protein',
                                'five_prime_UTR',
                                'three_prime_UTR',
                                'pcr_product',
                                'CDS',
                                'exon',
                                'intron'],
                        'GTF':['gene',
                               'mRNA',
                               'CDS',
                               'exon',
                               'start_codon',
                               'stop_codon']}
    
    return GFF_parent_check_level_1,GFF_parent_check_level_2,GTF_parent_check_level_1,GTF_parent_check_level_2,expected_feature_counts,expected_features

(
GFF_parent_check_level_1,
GFF_parent_check_level_2,
GTF_parent_check_level_1,
GTF_parent_check_level_2,
expected_feature_counts,
expected_features,
) = EXPECTED_DATA()

def test_fasta2oneline():
    GFFutils.fasta_seqs_to_oneline('multiline.fa','multiline.fa.tmp')
    observed = open('multiline.fa.tmp').read()

    expected = open('singleline.fa').read()
    print expected
    print observed
    assert observed == expected
    os.unlink('multiline.fa.tmp')

def test_clean_gff():
    # test the "full" cleaning -- remove some featuretypes, do sanity-checking,
    # add chr
    GFFutils.clean_gff('dirty.gff',newfn='cleaned.tmp',featuretypes_to_remove=['pcr_product','protein'],addchr=True)
    observed = open('cleaned.tmp').readlines()
    expected = open('fully-cleaned.gff').readlines()
    assert observed==expected
    os.unlink('cleaned.tmp')
    GFFutils.clean_gff('dirty.gff',featuretypes_to_remove=None, sanity_check=False)
    observed = open('dirty.gff.cleaned').read()
    expected = open('basic-cleaned.gff').read()
    assert observed == expected
    os.unlink('dirty.gff.cleaned')

def test_inspect_featuretypes():
    observed = GFFutils.inspect_featuretypes('FBgn0031208.gff')
    observed.sort()
    expected = ['CDS', 'exon', 'five_prime_UTR', 'gene', 'intron', 'mRNA', 'pcr_product', 'protein', 'three_prime_UTR']
    print observed
    print expected
    assert observed == expected

class GenericFeature(object):
    def setup(self):
        self.featureclass = self.__class__.featureclass
        if self.featureclass == 'GFF':
            self.Feature = getattr(GFFutils,'GFFFeature')
        if self.featureclass == 'GTF':
            self.Feature = getattr(GFFutils,'GTFFeature')

    def test_feature_creation(self):
        "Creation of GFF/GTFFeature from kwargs"

        # When strvals is False, if start and stop are None then raise
        # TypeError (cause you can't convert to int from None)
        kwargs = dict(chrom=None,
                      start=None,
                      stop=None,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name=None,
                      strvals=False)
        nt.assert_raises(TypeError,self.Feature,**kwargs)
       

        # But it's OK if they're None if strvals is True, since None -> 'None'
        kwargs = dict(chrom=None,
                      start=None,
                      stop=None,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name=None,
                      strvals=True)
        feature = self.Feature(**kwargs)
       
        # Test chrom retrieval
        kwargs = dict(chrom='chr2L',
                      start=1,
                      stop=10,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name=None,
                      strvals=True)
        feature = self.Feature(**kwargs)
        assert feature.chr == 'chr2L'
        assert feature.chrom == 'chr2L'
        assert feature.name == 'chr2L'
        
        # Test "chrom" retrieval if you specify "name" kwarg
        kwargs = dict(chrom=None,
                      start=1,
                      stop=10,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name='chr2L',
                      strvals=True)
        feature = self.Feature(**kwargs)
        assert feature.chr == 'chr2L'
        assert feature.chrom == 'chr2L'
        assert feature.name == 'chr2L'
        
        # specifying chrom and name should raise ValueError
        kwargs = dict(chrom='chr2L',
                      start=1,
                      stop=10,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name='chr2L',
                      strvals=False)
        nt.assert_raises(ValueError,self.Feature,**kwargs)

        kwargs = dict(chrom='chr2L',
                      start=1,
                      stop=10,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      name=None,
                      strvals=False)

        feature = self.Feature(**kwargs)
        assert len(feature.attributes._attrs)==0

        # should auto-wrap values into strings.
        feature.add_attribute('length',100)
        assert feature.attributes.length[0] == 100
        print 'feature.attributes._attrs:',feature.attributes._attrs
        print 'feature._strattributes:', feature._strattributes

        if self.featureclass == 'GFF':
            assert feature._strattributes == 'length=100;'
        if self.featureclass == 'GTF':
            assert feature._strattributes == 'length "100";'

        feature.add_attribute('ID','number1!')
        print 'current feature._strattributes:', feature._strattributes
        assert feature.attributes.ID[0] =='number1!'
        assert feature.attributes.length[0] == 100

        print 'current feature._strattributes:', feature._strattributes
        if self.featureclass == 'GFF':
            assert feature._strattributes == 'length=100;ID=number1!;'
        if self.featureclass == 'GTF':
            assert feature._strattributes == 'length "100";ID "number1!";'
        
        feature.remove_attribute('ID')
        print 'current feature._strattributes:', feature._strattributes
        if self.featureclass == 'GFF':
            assert feature._strattributes == 'length=100;'
        if self.featureclass == 'GTF':
            assert feature._strattributes == 'length "100";'

        feature.remove_attribute('length')
        print 'current feature._strattributes:', feature._strattributes
        assert feature._strattributes == ''

    def test_to_bed(self):
        kwargs = dict(chrom=None,
                      start=None,
                      stop=None,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      strvals=True)
        feature = self.Feature(**kwargs)
        nt.assert_raises(ValueError,feature.to_bed)
        feature.chrom = 'chr2L'
        feature.start = 10
        feature.stop = 20
        feature.strand = '+'
        assert feature.to_bed() == 'chr2L\t9\t20\n'
        assert feature.to_bed(6) == 'chr2L\t9\t20\tNone\tNone\t+\n'

        feature.value = 0
        feature.add_attribute('ID','new-feature')
        print feature.to_bed(6)
        assert feature.to_bed(6) == 'chr2L\t9\t20\tnew-feature\t0\t+\n'
    
    def test_length(self):
        kwargs = dict(chrom='chr3',
                      start=1,
                      stop=30,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      strvals=True)
        feature = self.Feature(**kwargs)
        assert len(feature) == 29

        # built-in len will complain if result is negative
        feature.start = 31
        nt.assert_raises(ValueError,len,feature)

    def test_TSS(self):
        kwargs = dict(chrom='chr3',
                      start=1,
                      stop=30,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      strvals=True)
        feature = self.Feature(**kwargs)

        # not stranded, so TSS should give ValueError
        args = (feature, 'TSS')
        nt.assert_raises(ValueError, getattr, *args)

        # Add a strand, and it should work.
        feature.strand = '+'
        assert feature.TSS == 1
        
        # Minus strand, too
        feature.strand = '-'
        assert feature.TSS == 30

    def test_midpoint(self):
        kwargs = dict(chrom='chr3',
                      start=1,
                      stop=30,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=None,
                      strvals=True)
        feature = self.Feature(**kwargs)
        print feature.midpoint
        assert feature.midpoint == 15

        feature.start = None
        args = (feature, 'midpoint')
        nt.assert_raises(ValueError, getattr, *args)

    def test_tostring(self):
        kwargs = dict(chrom='chr3',
                      start=1,
                      stop=30,
                      value=99,
                      strand='+',
                      phase=None,
                      attributes=None,
                      strvals=False)
        feature = self.Feature(**kwargs)
        print feature.tostring()
        assert feature.tostring() == 'chr3\t.\t.\t1\t30\t99.0\t+\t.\n'

        if self.featureclass == 'GFF':
            attributes = 'asdf=100'
        if self.featureclass == 'GTF':
            attributes = 'asdf "100"'
        kwargs = dict(chrom='chr2L',
                      start=1,
                      stop=10,
                      value=None,
                      strand=None,
                      phase=None,
                      attributes=attributes,
                      strvals=True)
        feature = self.Feature(**kwargs)

        # TODO: decide whether you want things converted to float or not in the attributes....
        if self.featureclass == 'GFF':
            expected = 'chr2L\t.\t.\t1\t10\t.\t.\t.\tasdf=100;\n'
        if self.featureclass == 'GTF':
            # no float conversion (not implemented)
            expected = 'chr2L\t.\t.\t1\t10\t.\t.\t.\tasdf "100";\n'
            
            # float converstion (currently implemented)
            expected = 'chr2L\t.\t.\t1\t10\t.\t.\t.\tasdf "100.0";\n'
        print 'observed:', feature.tostring()
        print 'expected:', expected
        assert feature.tostring() == expected

    def test_eq(self):
        kwargs = dict(chrom='chr2L',
                      start=1,
                      stop=10,
                      value=0.0,
                      strand='+',
                      featuretype='gene',
                      phase=None,
                      attributes=None,
                      strvals=False)
        original = self.Feature(**kwargs)

        def substitute(**newkwargs):
            kwargs_copy = copy.deepcopy(kwargs)
            for key,val in newkwargs.items():
                kwargs_copy[key]=val
                print key,val
            other = self.Feature(**kwargs_copy)
            return other
        
        # each call to substitute will substitute the kwargs to Feature
        # creation.

        # test the substitute() function
        other = substitute(chrom='chr2f')
        assert original != other

        
        other = substitute(chrom='chr2L')
        assert original == other


        assert original != substitute(chrom='fake')
        assert original != substitute(start=5)
        assert original != substitute(value=100)
        assert original != substitute(strand='-')

        if self.featureclass == 'GTF':
            assert original != substitute(attributes='asdf "one";')
        if self.featureclass == 'GFF':
            assert original != substitute(attributes='asdf=one;')

class TestGFFFeature(GenericFeature):
    featureclass = 'GFF'

class TestGTFFeature(GenericFeature):
    featureclass = 'GTF'

class GenericDBClass(object):
    featureclass = None
    def setup(self):
        """
        Creates a new GFFDB or GTFDB (depending on self.__class__.featureclass)
        """
        if os.path.exists(testdbfn):
            os.remove(testdbfn)

        self.featureclass = self.__class__.featureclass

        if self.featureclass == 'GFF':
            extension = '.gff'
            self.Feature = GFFutils.GFFFeature
        if self.featureclass == 'GTF':
            extension = '.gtf'
            self.Feature = GFFutils.GTFFeature
            
        self.fn = 'FBgn0031208' + extension
        self.dbfn = testdbfn
        if self.featureclass == 'GFF':
            #GFFutils.create_gffdb(self.fn, self.dbfn, verbose=False)
            #self.G = GFFutils.GFFDB(self.dbfn)
            creation_class = GFFutils.GFFDBCreator(self.dbfn)
            creation_class.create_from_file(self.fn)
            self.G = GFFutils.GFFDB(self.dbfn)
        if self.featureclass == 'GTF':
            #GFFutils.create_gtfdb(self.fn, self.dbfn, verbose=False)
            creation_class = GFFutils.GTFDBCreator(self.dbfn)
            creation_class.create_from_file(self.fn)
            self.G = GFFutils.GTFDB(self.dbfn)

        self.conn = sqlite3.connect(self.dbfn)
        self.c = self.conn.cursor()

    def table_test(self):
        """Right tables exist?"""
        expected_tables = ['features','relations']
        expected_tables.sort()
        self.c.execute('select name from sqlite_master where type="table"')
        observed_tables = [str(i[0]) for i in self.c]
        observed_tables.sort()
        print expected_tables
        print observed_tables
        assert expected_tables == observed_tables
    
    def _count1(self,featuretype):
        """
        Count using SQL
        """
        self.c.execute('select count() from features where featuretype = ?',(featuretype,))
        results = self.c.fetchone()[0]
        print 'count1("%s") says: %s' % (featuretype,results)
        return results

    def _count2(self,featuretype):
        """
        Count GFF lines.
        """
        cnt = 0
        for line in open(self.fn):
            if line.startswith('#'):
                continue
            L = line.split()
         
            if len(L) < 3:
                continue
     
            if L[2] == featuretype:
                cnt += 1
        print 'count2("%s") says: %s' % (featuretype, cnt)
        return cnt
    
    def _count3(self,featuretype):
        """
        Count with the count_features_of_type method.
        """
        results = self.G.count_features_of_type(featuretype)
        print 'count3("%s") says: %s' % (featuretype, results)
        return results
    
    def _count4(self,featuretype):
        """
        Count by iterating over all features of this type
        """
        cnt = 0
        for i in self.G.features_of_type(featuretype):
            cnt += 1
        print 'count4("%s") says: %s' % (featuretype,cnt)
        return cnt
    
    def featurecount_test(self):
        """
        Right number of each featuretype, using 3 different ways of
        counting?
        """
        featuretypes = ['gene',
                        'mRNA',
                        'CDS',
                        'exon',
                        'five_prime_UTR',
                        'three_prime_UTR',
                        'pcr_product',
                        'protein','intron']
        for featuretype in featuretypes:
            rawsql_cnt = self._count1(featuretype)
            count_feature_of_type_cnt = self._count3(featuretype)
            iterator_cnt = self._count4(featuretype)
            try:
                hard_count = expected_feature_counts[self.featureclass][featuretype]
            except KeyError:
                hard_count = 0

            # count2 is not an appropriate test for GTF files, since "gene" is
            # not explicitly listed -- only implied by the boundaries of CDSs:
            #
            if self.featureclass == 'GFF':
                gffparsed_cnt = self._count2(featuretype)
                assert rawsql_cnt == count_feature_of_type_cnt == iterator_cnt == gffparsed_cnt == hard_count
            
            if self.featureclass == 'GTF':
                print 'hard count:',hard_count 
                assert rawsql_cnt == count_feature_of_type_cnt == iterator_cnt == hard_count

    def total_features_test(self):
        """did all features in the GFF file make it into the database?"""

        # Doesn't make sense for GTF, which has gene features implied rather
        # than explicitly stated.
        if self.featureclass == 'GTF':
            return
        cnt = 0
        if self.__class__.featureclass == 'GFF':
            iterator = GFFutils.GFFFile(self.fn)
        if self.__class__.featureclass == 'GTF':
            iterator = GFFutils.GTFFile(self.fn)
        for feature in iterator:
            cnt += 1
        self.c.execute('select count() from features')
        total = self.c.fetchone()[0]
        assert cnt == total

    def parents_test(self):
        "DB retrieval of parents matches expected?"
        # Checks the hand-entered lookup dicts at the beginning of this module
        # to make sure parents match up.
        if self.featureclass == 'GFF':
            parents1 = GFF_parent_check_level_1
            parents2 = GFF_parent_check_level_2
        if self.featureclass == 'GTF':
            parents1 = GTF_parent_check_level_1
            parents2 = GTF_parent_check_level_2
        for child, expected_parents in parents1.items():
            observed_parents = [i.id for i in self.G.parents(child, level=1)]
            print 'observed parents for %s:' % child, set(observed_parents)
            print 'expected parents for %s:' % child, set(expected_parents)
            assert set(observed_parents) == set(expected_parents)

    def identity_test(self):
        """
        Actual string in file identical to tostring() method output?
        """

        for i in open(self.fn):
            if i.startswith('#'):
                continue
            first_line = i
            break
        if self.featureclass == 'GFF':
            iterator = iter(GFFutils.GFFFile(self.fn))
        if self.featureclass == 'GTF':
            iterator = iter(GFFutils.GTFFile(self.fn))
            
        first_feature = iterator.next()
        feature_string = first_feature.tostring()

        # We're always returning the attributes with a trailing ";"
        if first_line.endswith(';\n'):
            print 'first line:',first_line
            print 'first feat:',first_feature.tostring()
            assert feature_string == first_line
        else:
            first_line = first_line.strip()+';'+'\n'
            print 'first line:',first_line
            print 'first feat:',first_feature.tostring()
            assert feature_string == first_line

    def all_test(self):
        """all() returns expected number of features?"""
        results = list(self.G.all())
        if self.featureclass == 'GFF':
            print len(results)
            assert len(results) == 27
        if self.featureclass == 'GTF':
            print len(results)

            # 23 entries;
            #  3 genes, but 1 is noncoding
            #  4 transcripts (2 for one gene, one for the fake coding and one for the fake non-coding 
            #  6 exons (2 in one isoform, 3 in the other (one of which is actually shared), 1 for each fake)
            #  5 CDS (cause one gene is non-coding)
            #  2 start_codons (1 shared between isoforms, one for fake coding)
            #  3 stop_codons (1 for each coding transcript)
            assert len(results) == 23  

    def features_of_type_test(self):
        """features_of_type() returns expected number of features?"""
        d = expected_feature_counts[self.featureclass]
        for key,val in d.items():
            observed = len(list(self.G.features_of_type(key)))
            print 'key:',"'"+key+"'",'val:',val,'observed:', observed
            assert observed == val
        
        # now either catch all...
        for key,val in d.items():
            observed = len(list(self.G.features_of_type(key, chrom='chr2L', start=1, stop=100000)))
            assert observed == val
        
        # or none...
        for key,val in d.items():
            
            # wrong chrom should return 0
            observed = len(list(self.G.features_of_type(key, chrom='chrX', strand='+', start=1, stop=10000)))
            assert observed == 0
            
            # wrong strand should return 0
            observed = list(self.G.features_of_type(key, chrom='chr2L', strand='-', start=9999, stop=100000))
            expected_len = 2
            if key == 'gene' or key == 'mRNA':
                print 'observed:',observed
                print 'expected len:',expected_len
                assert len(observed) == expected_len

            
            # too far into chrom should return 0
            observed = len(list(self.G.features_of_type(key, chrom='chr2L', strand='+', start=100000, stop=1e6)))
            assert observed == 0
            
            # reversed start/stop should return 0
            observed = len(list(self.G.features_of_type(key, chrom='chr2L', start=10000, stop=1)))
            assert observed == 0

    def features_test(self):
        """All featuretypes present and accounted for?"""
        expected = expected_features[self.featureclass]
        results = list(self.G.features())
        assert set(results) == set(expected)

    def strand_test(self):
        """Expected strands in db?"""
        assert set(self.G.strands()) == set(['-','+'])

    def chrom_test(self):
        """Expected chromosomes in db?"""
        # Get the chromosomes that are in the GFF/GTF file, and make sure
        # they made it into the database.
        gffchroms = []
        for line in open(self.fn):
            if line.startswith('>'):
                break
            if line.startswith('#'):
                continue
            L = line.split()
            if len(L) < 1:
                continue
            gffchroms.append(L[0])
        print 'observed:',set(self.G.chromosomes())
        print 'expected:',set(gffchroms)
        assert set(gffchroms) == set(self.G.chromosomes())

    def closest_features_test(self):
        """Expected closest features returned?"""

        # closest one to the beginning, on plus strand
        observed_dist, observed_feature = self.G.closest_feature(chrom='chr2L', 
                                                                 pos=1,
                                                                 featuretype='gene',
                                                                 strand='+')
        print observed_feature
        assert observed_feature.id == 'FBgn0031208'
        
        # closest one to beginning, ignoring the first one.  Should return None.
        result = self.G.closest_feature(chrom='chr2L', 
                                        pos=1, 
                                        featuretype='gene',
                                        strand='+', 
                                        ignore=['FBgn0031208'])
        print result
        assert result == (None,None)
        
        # closest mRNA to beginning, on - strand.
        observed_dist, observed_feature = self.G.closest_feature(chrom='chr2L', 
                                                                 pos=1,
                                                                 featuretype='mRNA',
                                                                 strand='-', 
                                                                 ignore=None,)
        print observed_feature
        assert observed_feature.id == 'transcript_Fk_gene_1'

        # TODO: the GFF/GTF files probably need some more genes added to be
        # able to test the "upstream" # and "downstream" functionality.
        
    
    def feature_not_found_test(self):
        """Correct exception thrown upon asking for nonexistent feature?"""
        try:
            self.G['i am not a real feature']
        except GFFutils.FeatureNotFoundError:
            pass
        
        try:
            self.G['another fake one']
        except GFFutils.FeatureNotFoundError as e:
            print e
            assert str(e) == 'another fake one'

    def getitem_and_eq_test(self):
        # test access of a string or a feature, make sure you get back the same
        # thing.
        if self.featureclass == 'GFF':
            expected = 'chr2L\tFlyBase\tgene\t7529\t9484\t.\t+\t.\tID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0008234,GO:0006508;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB:AE003590,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q6KEV3,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,BIOGRID:59420,dedb:5870,GenomeRNAi_gene:33155,ihop:59383;derived_computed_cyto=21A5-21A5;gbunit=AE014134;\n'
        if self.featureclass == 'GTF':
            expected = """chr2L	protein_coding	exon	7529	8116	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	exon	7529	8116	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	start_codon	7680	7682	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	start_codon	7680	7682	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	CDS	7680	8116	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	CDS	7680	8116	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	CDS	8193	8589	.	+	2	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	exon	8193	8589	.	+	2	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	CDS	8193	8610	.	+	2	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	exon	8193	9484	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	stop_codon	8611	8613	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
chr2L	protein_coding	CDS	8668	9276	.	+	2	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	exon	8668	9484	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
chr2L	protein_coding	stop_codon	9277	9279	.	+	0	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
"""
        observed_feature_1 = self.G['FBgn0031208']

        if self.featureclass == 'GTF':
            nt.assert_raises(ValueError,observed_feature_1.tostring)
            observed_string_1 = observed_feature_1.tostring(gtfdb=self.G,include_id=False)

        if self.featureclass == 'GFF':
            observed_string_1 = observed_feature_1.tostring()

        #print 'expected :',expected
        #print 'observed1:',observed_string_1

        #d = difflib.ndiff(expected.splitlines(True), observed_string_1.splitlines(True))
        #print ''.join(d)
        print observed_string_1
        assert observed_string_1 == expected
        
        observed_feature_2 = self.G[observed_feature_1]
        
        if self.featureclass == 'GFF':
            observed_string_2 = observed_feature_2.tostring()
        if self.featureclass == 'GTF':
            observed_string_2 = observed_feature_2.tostring(gtfdb=self.G, include_id=False)
        print 'expected :',expected
        print 'observed2:',observed_string_2
        assert observed_string_2 == expected

        assert observed_feature_1 == observed_feature_2

    def overlapping_features_test(self):
        # test everything against a hard-coded list.
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L', start=7529, stop=9484,completely_within=False)])
        print observed
        if self.featureclass == 'GFF':
            expected_hardcoded = ['CDS_FBgn0031208:1_737',
                             'CDS_FBgn0031208:2_737',
                             'CDS_FBgn0031208:3_737',
                             'CDS_FBgn0031208:4_737',
                             'FBgn0031208',
                             'FBgn0031208:2',
                             'FBgn0031208:3',
                             'FBgn0031208:4',
                             'FBpp0289913',
                             'FBpp0289914',
                             'FBtr0300689',
                             'FBtr0300690',
                             'INC121G01_pcr_product',
                             'exon:chr2L:7529-8116',
                             'five_prime_UTR_FBgn0031208:1_737', 
                             'intron_FBgn0031208:1_FBgn0031208:2',
                             'intron_FBgn0031208:1_FBgn0031208:3', 
                             'intron_FBgn0031208:2_FBgn0031208:4',
                             'three_prime_UTR_FBgn0031208:3_737',
                             'three_prime_UTR_FBgn0031208:4_737', 
                             ]
        if self.featureclass == 'GTF':
            expected_hardcoded = ['CDS:chr2L:7680-8116',
                                  'CDS:chr2L:8193-8610',
                                  'CDS:chr2L:8193-8589', 
                                  'CDS:chr2L:8668-9276',
                                  'FBgn0031208',
                                  'FBtr0300689',
                                  'FBtr0300690', 
                                  'exon:chr2L:7529-8116',
                                  'exon:chr2L:8193-9484',
                                  'exon:chr2L:8193-8589',
                                  'exon:chr2L:8668-9484', 
                                  'start_codon:chr2L:7680-7682', 
                                  'stop_codon:chr2L:8611-8613',
                                  'stop_codon:chr2L:9277-9279']

        print 'observed:',set(observed)
        print 'expected:',set(expected_hardcoded)
        assert set(observed) == set(expected_hardcoded)

        # assert that features that overlap the gene contain the gene's children and grandchildren
        observed = [i.id for i in self.G.overlapping_features(chrom='chr2L', start=7529, stop=9484, featuretype='gene')]
        observed += [i.id for i in self.G.overlapping_features(chrom='chr2L', start=7529, stop=9484, featuretype='exon')]
        observed += [i.id for i in self.G.overlapping_features(chrom='chr2L', start=7529, stop=9484, featuretype='intron')]
        observed += [i.id for i in self.G.overlapping_features(chrom='chr2L', start=7529, stop=9484, featuretype='mRNA')]
        observed.sort()

        expected_generated = list(self.G.children('FBgn0031208',level=1,featuretype='mRNA'))
        expected_generated += list(self.G.children('FBgn0031208',level=2,featuretype='exon'))
        expected_generated += list(self.G.children('FBgn0031208',level=2,featuretype='intron'))
        expected_generated.append(self.G['FBgn0031208'])
        expected_generated = sorted([i.id for i in expected_generated])
        print observed
        print expected_generated
        print set(observed).difference(expected_generated)
        assert observed == expected_generated
        
        # check that the first exon, which is chr2L:7529-8116, overlaps lots of stuff....
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L',start=7529,stop=8116, completely_within=False, strand='+')])
        print 'observed:',observed
        if self.featureclass == 'GFF':
            expected = ['CDS_FBgn0031208:1_737',
                        'FBgn0031208', 
                        'FBpp0289913', 
                        'FBpp0289914',
                        'FBtr0300689',
                        'FBtr0300690', 
                        'exon:chr2L:7529-8116',
                        'five_prime_UTR_FBgn0031208:1_737', 
                        ]
        if self.featureclass == 'GTF':
            # things like the last exons and stop codons shouldn't be overlapping.
            expected = ['CDS:chr2L:7680-8116', 
                        'FBgn0031208',
                        'FBtr0300689',
                        'FBtr0300690',
                        'exon:chr2L:7529-8116',
                        'start_codon:chr2L:7680-7682']
        print 'expected:',expected
        assert observed == expected

        # on the opposite strand, though, you should have nothing.
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L',start=7529,stop=8116, completely_within=False, strand='-')])
        assert observed == []


        # completely_within, however, should only return the UTR, exon, and CDS for GFF, (and exon, CDS, and start codon for GTF)
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L',start=7529,stop=8116, completely_within=True)])
        print 'observed:',observed
        if self.featureclass == 'GFF':
            expected = ['CDS_FBgn0031208:1_737','exon:chr2L:7529-8116','five_prime_UTR_FBgn0031208:1_737',]
        if self.featureclass == 'GTF':
            # multiple exons, multiple CDSs, multiple start_codons.  Still not
            # sure if this is the best way to do it, since it's not identical
            # to GFF output.
            expected = ['CDS:chr2L:7680-8116',
                        'exon:chr2L:7529-8116',
                        'start_codon:chr2L:7680-7682']
        print 'expected:',expected
        assert observed == expected

        # stop truncated by 1 bp should only get you the UTR for GFF, or just the start codon[s] for GTF.
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L',start=7529,stop=8115, completely_within=True)])
        print 'observed:',observed
        if self.featureclass == 'GFF':
            expected = ['five_prime_UTR_FBgn0031208:1_737']
        if self.featureclass == 'GTF':
            expected = ['start_codon:chr2L:7680-7682']
        print 'expected:',expected
        assert observed == expected
        
        # and the 5'UTR, but with start truncated 1 bp, should get you nothing -- for both GFF and GTF.
        observed = sorted([i.id for i in self.G.overlapping_features(chrom='chr2L',start=7530,stop=7679, completely_within=True)])
        print 'observed:',observed
        expected = []
        print 'expected:',expected
        assert observed == expected
    
    def merge_features_test(self):
        features_to_merge = []
        chrom = 'chr2L'
        coords = [(1,10,'exon'), (5,50,'intron'), (100,105,'CDS')]
        for start, stop, featuretype in coords:
            feature = self.Feature(chrom=chrom,start=start,stop=stop,featuretype=featuretype)
            features_to_merge.append(feature)
        observed = list(self.G.merge_features(features_to_merge))
        expected = [self.Feature(chrom=chrom,start=1,stop=50,featuretype='merged_CDS_exon_intron',strand=None),
                    self.Feature(chrom=chrom,start=100,stop=105,featuretype='merged_CDS_exon_intron',strand=None)]

        print 'observed:', observed
        print 'expected:',expected 
        assert observed == expected

        # make sure you complain if more than one strand and ignore_strand=False
        diff_strands = [self.Feature(chrom=chrom,start=1,stop=10,strand='+'),
                        self.Feature(chrom=chrom,start=5,stop=20,strand='-')
                       ]
        merged = self.G.merge_features(diff_strands)
        nt.assert_raises(ValueError, list, merged)
        
        # But this should run without raising an error
        merged = list(self.G.merge_features(diff_strands,ignore_strand=True))
        
        diff_chroms = [self.Feature(chrom='chr2L',start=1,stop=10,strand='+'),
                       self.Feature(chrom='chr3R',start=5,stop=20,strand='+')
                      ]
        merged = self.G.merge_features(diff_chroms)
        nt.assert_raises(NotImplementedError, list, merged)
        exons = self.G.children('FBgn0031208',featuretype='exon', level=2)
        merged_exons = self.G.merge_features(exons)
        observed = list(merged_exons)
        
        expected = [
                    self.Feature(chrom='chr2L',start=7529,stop=8116, strand='+', featuretype='merged_exon'),
                    self.Feature(chrom='chr2L',start=8193,stop=9484, strand='+', featuretype='merged_exon'),
                   ]

        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

    def interfeatures_test(self):
        chrom = 'chr2L'
        coords = [(1,10,'exon'), (50,100,'exon'), (200,300,'other')]
        features = []
        for start,stop,featuretype in coords:
            features.append(self.Feature(chrom=chrom,start=start,stop=stop,featuretype=featuretype))
        observed = list(self.G.interfeatures(features))
        
        expected = [self.Feature(chrom=chrom,start=11,stop=49,featuretype='inter_exon_exon'),
                    self.Feature(chrom=chrom,start=101,stop=199,featuretype='inter_exon_other'),]
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected
    
    def execute_test(self):
        query = 'SELECT id, chrom, start, stop FROM features WHERE id = "FBgn0031208"'
        expected = [('FBgn0031208','chr2L',7529,9484)]
        observed = list(self.G.execute(query))
        print observed
        print expected
        assert observed == expected
            

    def reflat_test(self):
        gene_id = 'FBgn0031208'
        gene = self.G[gene_id]
        observed = ''.join(self.G.refFlat(gene))
        if self.featureclass == 'GFF':
            expected = """CG11023-RB	FBgn0031208	chr2L	+	7529	9484	7680	8610	2	7529,8193,	8116,9484,
CG11023-RC	FBgn0031208	chr2L	+	7529	9484	7680	9276	3	7529,8193,8668,	8116,8589,9484,
"""
        if self.featureclass == 'GTF':
            expected = """FBtr0300689	FBgn0031208	chr2L	+	7529	9484	7680	8610	2	7529,8193,	8116,9484,
FBtr0300690	FBgn0031208	chr2L	+	7529	9484	7680	9276	3	7529,8193,8668,	8116,8589,9484,
"""
        print 'observed:'
        print observed
        print 'expected:'
        print expected
        assert observed == expected

    def promoter_test(self):
        # test defaults of bidirectional=True and dist=1000
        observed = self.G.promoter(id='FBgn0031208')
        expected = self.Feature(chrom='chr2L',start=6529,stop=8529,strand='+')
        expected.add_attribute('ID','promoter:FBgn0031208')
        expected.featuretype = 'promoter'
        print 'observed:'
        print observed
        print 'len(observed):',len(observed)
        print 'expected:'
        print expected
        assert len(observed) == 2000
        assert observed == expected


        # Test various kwargs
        observed = self.G.promoter(id='FBgn0031208',dist=100,bidirectional=False)
        expected = self.Feature(chrom='chr2L',start=7429,stop=7529,strand='+')
        expected.add_attribute('ID','promoter:FBgn0031208')
        expected.featuretype = 'promoter'
        print 'observed:'
        print observed
        print 'len(observed):',len(observed)
        print 'expected:'
        print expected
        assert len(observed) == 100
        assert observed == expected

        # Test truncation
        # TODO: add another gene in GTF file....
        if self.featureclass == 'GFF':
            observed = self.G.promoter(id='Fk_gene_1',dist=5000,bidirectional=False,truncate_at_next_feature='gene')
            expected = self.Feature(chrom='chr2L',start=11000,stop=11500,strand='-')
            expected.add_attribute('ID','promoter:Fk_gene_1')
            expected.featuretype = 'promoter'
            print 'observed:'
            print observed
            print 'len(observed):',len(observed)
            print 'expected:'
            print expected
            assert len(observed) == 500
            assert observed == expected
    
    def attribute_search_test(self):
        if self.featureclass == 'GFF':
            observed = self.G.attribute_search('FBan0011023')
            expected = self.G['FBgn0031208']
            assert len(observed)==1
            assert observed[0] == expected

        if self.featureclass == 'GTF':
            observed = self.G.attribute_search('Fk_gene_1')
            expected = self.G['Fk_gene_1']
            assert len(observed)==1
            assert observed[0] == expected

    def exonic_bp_test(self):
        observed = self.G.exonic_bp('FBgn0031208')
        expected = 1878
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

        # If no exons, then exonic length is zero
        observed = self.G.exonic_bp('transcript_Fk_gene_1')
        expected = 0
        assert observed == expected
    
    def random_feature_test(self):
        # really all we can do is make sure the featuretype is correct and that
        # upon retrieving the same feature from the DB, we get identity.
        random_gene = self.G.random_feature('gene')
        assert random_gene.featuretype == 'gene'
        assert self.G[random_gene.id] == random_gene
    
    def coding_genes_test(self):
        observed = sorted([i.id for i in self.G.coding_genes()])
        expected = ['FBgn0031208','Fk_gene_1']
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected
    
    def n_exon_isoforms_test(self):
        exons = list(self.G.children('FBgn0031208',level=2, featuretype='exon'))
        exons.sort(key=lambda x: x.start)
        exon = exons[0]
        observed = self.G.n_exon_isoforms(exon)
        expected = 2
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

    def exons_gene_test(self):
        # exons_gene uses a different query than the [already-tested] children.
        # So hopefully this isn't a tautology
        exons = list(self.G.children('FBgn0031208',level=2, featuretype='exon'))
        exons.sort(key=lambda x: x.start)
        exon = exons[0]
        observed = self.G.exons_gene(exon)
        expected = 'FBgn0031208'
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected
    
    def teardown(self):
        pass
        #os.unlink(testdbfn)

class TestGFFDBClass(GenericDBClass):
    featureclass = 'GFF'

class TestGTFDBClass(GenericDBClass):
    featureclass = 'GTF'
    def UTR_(self):
        
        observed = self.G.UTRs(self.G['FBtr0300689'])
        expected = [
                    self.Feature(chrom='chr2L',start=7529,stop=7679,strand='+',featuretype='five_prime_UTR'),
                    self.Feature(chrom='chr2L',start=8611,stop=9484,strand='+',featuretype='three_prime_UTR'),
                   ]
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

        observed = self.G.UTRs(self.G['FBtr0300690'])
        expected = [
                    self.Feature(chrom='chr2L',start=7529,stop=7679,strand='+',featuretype='five_prime_UTR'),
                    self.Feature(chrom='chr2L',start=9277,stop=9484,strand='+',featuretype='three_prime_UTR'),
                   ]
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

class TestGenomeWithGFFDB(object):
    def setup(self):
        if os.path.exists('seqtest.db'):
            os.unlink('seqtest.db')
        self.genome = GFFutils.Genome('singleline.fa',debug=True)
        GFFutils.create_gffdb('seqtest.gff','seqtest.db')
        self.gffdb = GFFutils.GFFDB('seqtest.db')
    def teardown(self):
        if os.path.exists('seqtest.db'):
            os.unlink('seqtest.db')
    def test_spliced_transcript(self):
        featureid = 'transcript1'
        observed = self.genome.spliced_transcript(featureid=featureid,
                                                   gffdb=self.gffdb,
                                                   featuretype='exon')
        expected = 'GAAGCAGAACAGATATTTAGATTGCCTCTCATTTATAGGGAGAAATATGATCGCGTATGCGAGTGCCAACATATTGTGCTCTTTGATTTTT'
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

        feature = self.gffdb[featureid]
        observed = self.genome.spliced_transcript(featureid=feature,
                                                   gffdb=self.gffdb,
                                                   featuretype='exon')
        expected = 'GAAGCAGAACAGATATTTAGATTGCCTCTCATTTATAGGGAGAAATATGATCGCGTATGCGAGTGCCAACATATTGTGCTCTTTGATTTTT'
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected
    
    def test_splice_junctions(self):
        featureid = 'transcript1'
        observed = self.genome.splice_junctions(featureid=featureid,
                                                gffdb=self.gffdb,
                                                padding=5)
        observed = list(observed)
        expected = ['CTCATTTATAGG','TATGCGAGTGCC']
        print 'observed',observed
        print 'expected',expected
        assert observed == expected
        
        # Make sure it works for the feature, too.
        feature = self.gffdb[featureid]
        observed = self.genome.splice_junctions(featureid=feature,
                                                gffdb=self.gffdb,
                                                padding=5)
        observed = list(observed)
        print 'observed',observed
        print 'expected',expected
        assert observed == expected

class TestGenome(object):
    def setup(self):
        self.genome = GFFutils.Genome('singleline.fa',debug=True)

    def test_chrom_lens(self):
        observed = self.genome.chromlens['chrX']
        expected = 700
        print 'observed:', observed, 'expected:', expected
        assert observed == expected

        observed = self.genome.chromlens['chr2L']
        expected = 450
        print 'observed:', observed, 'expected:', expected
        assert observed == expected

    def test_chrom_pos(self):
        observed = self.genome.startinds['chrX']
        expected = 464
        print 'observed:', observed, 'expected:', expected
        assert observed == expected

        observed = self.genome.startinds['chr2L']
        expected = 7
        print 'observed:', observed, 'expected:', expected
        assert observed == expected

    def test_sequence(self):
        observed = self.genome.sequence('chr2L',1,10)
        expected = 'CGACAATGCA'
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected
        
        observed = self.genome.sequence('chrX',14,58)
        expected = 'ATGCCCACTGTGGGGAATTTACCAGCAGCCCGCACACTTAGCCGG'
        print 'observed:',observed
        print 'expected:',expected
        assert observed == expected

        # Check that non-existent chroms and out-of-range positions raise the
        # right kind of errors.
        args = ('chr99',1,10)
        nt.assert_raises(KeyError, self.genome.sequence, *args)
        
        args = ('chr2L',450,451)
        nt.assert_raises(ValueError, self.genome.sequence, *args)

        args = ('chrX',1000,2000)
        nt.assert_raises(ValueError, self.genome.sequence, *args)
        
        # Check the ends of the chroms
        observed = self.genome.sequence('chr2L',450,450)
        expected = 'A'
        print observed
        print expected
        assert observed == expected

        observed = self.genome.sequence('chrX',699,700)
        expected = 'AC'
        print observed 
        print expected
        assert observed == expected

    def test_revcomp_sequence(self):
        plus = self.genome.sequence('chr2L',36,43)
        minus = self.genome.sequence('chr2L',36,43,'-')
        print 'plus :',plus
        print 'minus:',minus
        assert plus == 'TAGATTGC'
        assert minus == 'GCAATCTA'

    
    def test_sequence_from_feature(self):
        gfffeature = GFFutils.GFFFeature(chrom='chr2L',
                                         start=36,
                                         stop=43,
                                         strand='+')
        observed = self.genome.sequence_from_feature(gfffeature)
        expected = 'TAGATTGC'
        assert observed == expected

        gfffeature = GFFutils.GFFFeature(chrom='chr2L',
                                         start=36,
                                         stop=43,
                                         strand='-')
        expected = 'GCAATCTA'
        observed = self.genome.sequence_from_feature(gfffeature)
        assert observed == expected       
    

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

