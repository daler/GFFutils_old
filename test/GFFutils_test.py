import GFFutils
import os
import sqlite3
import nose.tools as nt
import difflib
import pprint


testdbfn = 'testing.db'

def EXPECTED_DATA():
    # list the children and their expected first-order parents for the GFF test file.
    GFF_parent_check_level_1 = {'FBtr0300690':['FBgn0031208'],
                                'FBtr0300689':['FBgn0031208'],
                                'unnamed_exon_1':['FBtr0300689','FBtr0300690'],
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
                                'unnamed_exon_1':['FBgn0031208'],
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
                                'exon:FBtr0089256:7529-8116':['FBtr0089256'],
                                'FBtr0089256':['FBgn0031208'],
                                'CDS:FBtr0089256:7680-8116':['FBtr0089256'],
                                'start_codon:FBtr0089256:7680-7682':['FBtr0089256'],
                                'exon:FBtr0089256:8229-8589':['FBtr0089256'],
                                'CDS:FBtr0089256:8229-8589':['FBtr0089256'],
                                'exon:FBtr0089256:8668-9491':['FBtr0089256'],
                                'CDS:FBtr0089256:8668-9273':['FBtr0089256'],
                                'stop_codon:FBtr0089256:9274-9276':['FBtr0089256'],
                               }
    GTF_parent_check_level_2 = {
                                'exon:FBgn0031208:7529-8116':['FBgn0031208'],
                                'FBgn0031208':['FBgn0031208'],
                                'CDS:FBgn0031208:7680-8116':['FBgn0031208'],
                                'start_codon:FBgn0031208:7680-7682':['FBgn0031208'],
                                'exon:FBgn0031208:8229-8589':['FBgn0031208'],
                                'CDS:FBgn0031208:8229-8589':['FBgn0031208'],
                                'exon:FBgn0031208:8668-9491':['FBgn0031208'],
                                'CDS:FBgn0031208:8668-9273':['FBgn0031208'],
                                'stop_codon:FBgn0031208:9274-9276':['FBgn0031208'],
                               }

    expected_feature_counts = {
                'GFF':{'gene':2,
                       'mRNA':3,
                       'exon':5,
                       'CDS':5,
                       'five_prime_UTR':1,
                       'intron':3,
                       'pcr_product':1,
                       'protein':2,
                       'three_prime_UTR':2},
                'GTF':{'gene':2,
                       'mRNA':2,
                       'CDS':4,
                       'exon':4,
                       'start_codon':1,
                       'stop_codon':1}
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

def test_inspect_featuretypes():
    observed = GFFutils.inspect_featuretypes('FBgn0031208.gff')
    observed.sort()
    expected = ['CDS', 'exon', 'five_prime_UTR', 'gene', 'intron', 'mRNA', 'pcr_product', 'protein', 'three_prime_UTR']
    print observed
    print expected
    assert observed == expected

class TestGFFFeature(object):
    featureclass = 'GFF'

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
                      strvals=True)
        feature = self.Feature(**kwargs)
        assert feature.chr == 'chr2L'
        assert feature.chrom == 'chr2L'

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

class TestGTFFeature(TestGFFFeature):
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
        if self.featureclass == 'GTF':
            extension = '.gtf'
            
        self.fn = 'FBgn0031208' + extension
        self.dbfn = testdbfn
        if self.featureclass == 'GFF':
            GFFutils.create_gffdb(self.fn, self.dbfn)
            self.G = GFFutils.GFFDB(self.dbfn)
        if self.featureclass == 'GTF':
            GFFutils.create_gtfdb(self.fn, self.dbfn)
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
            assert len(results) == 24
        if self.featureclass == 'GTF':
            assert len(results) == 14
 
    def features_of_type_test(self):
        """features_of_type() returns expected number of features?"""
        d = expected_feature_counts[self.featureclass]
        for key,val in d.items():
            observed = len(list(self.G.features_of_type(key)))
            print key,val
            print observed
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
            observed = len(list(self.G.features_of_type(key, chrom='chr2L', strand='-', start=9999, stop=100000)))
            if key == 'gene' or key == 'mRNA':
                assert observed == 1

            
            # too far into chrom should return 0
            observed = len(list(self.G.features_of_type(key, chrom='chr2L', strand='+', start=10000, stop=1e6)))
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

    def teardown(self):
        os.unlink(testdbfn)

class TestGFFDBClass(GenericDBClass):
    featureclass = 'GFF'

class TestGTFDBClass(GenericDBClass):
    featureclass = 'GTF'

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
    
