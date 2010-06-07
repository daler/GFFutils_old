import GFFutils
import os
import sqlite3

testdbfn = 'testing.db'

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

class TestGFFDBClass(object):
    featureclass = 'GFF'
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
        if self.featureclass == 'GTF':
            first_feature.remove_attribute('ID')
        print 'first line:',first_line
        print 'first feat:',first_feature.tostring()
        feature_string = first_feature.tostring()
        if feature_string == first_line:
            return
        assert False

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
            if line.startswith('#'):
                continue
            L = line.split()
            if len(L) < 1:
                continue
            gffchroms.append(L[0])
        assert set(gffchroms) == set(self.G.chromosomes())

    def closest_features_test(self):
        """Expected closest features returned?"""

        # closest one to the beginning, on plus strand
        observed_dist, observed_id = self.G.closest_feature(chrom='chr2L', pos=1, featuretype='gene', strand='+', 
                                                            ignore=None, direction=None)
        assert observed_id == 'FBgn0031208'
        
        # closest one to beginning, ignoring the first one.  Should return None.
        result = self.G.closest_feature(chrom='chr2L', pos=1, featuretype='gene', strand='+', 
                                        ignore='FBgn0031208', direction=None)
        assert result is None
        
        # closest mRNA to beginning, on - strand.
        observed_dist, observed_id = self.G.closest_feature(chrom='chr2L', pos=1, featuretype='mRNA', strand='-', 
                                                            ignore=None, direction=None)

        assert observed_id == 'transcript_Fk_gene_1'

        # TODO: the GFF/GTF files probably need some more genes added to be
        # able to test the "upstream" # and "downstream" functionality.
        

    
    def teardown(self):
        pass
        os.remove(testdbfn)

class TestGTFDBClass(TestGFFDBClass):
    featureclass = 'GTF'

