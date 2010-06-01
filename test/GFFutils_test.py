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
                            'exon:FBtr0089256:1':['FBtr0089256'],
                            'FBtr0089256':['FBgn0031208'],
                            'CDS:FBtr0089256:1':['FBtr0089256'],
                            'start_codon:FBtr0089256:1':['FBtr0089256'],
                            'exon:FBtr0089256:2':['FBtr0089256'],
                            'CDS:FBtr0089256:2':['FBtr0089256'],
                            'exon:FBtr0089256:3':['FBtr0089256'],
                            'CDS:FBtr0089256:3':['FBtr0089256'],
                            'stop_codon:FBtr0089256:1':['FBtr0089256'],
                           }
GTF_parent_check_level_2 = {
                            'exon:FBtr0089256:3':['FBgn0031208'],
                            'CDS:FBtr0089256:3':['FBgn0031208'],
                            'stop_codon:FBtr0089256:1':['FBgn0031208'],
                            'start_codon:FBtr0089256:1':['FBgn0031208'],
                            'exon:FBtr0089256:1':['FBgn0031208'],
                            'CDS:FBtr0089256:1':['FBgn0031208'],
                            'exon:FBtr0089256:2':['FBgn0031208'],
                            'CDS:FBtr0089256:2':['FBgn0031208'],
                           }

class TestGFFDBClass(object):
    featureclass = 'GFF'
    def setup(self):
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
        """ensure the right tables exist"""
        expected_tables = ['features','relations']
        expected_tables.sort()
        self.c.execute('select name from sqlite_master where type="table"')
        observed_tables = [str(i[0]) for i in self.c]
        observed_tables.sort()
        print expected_tables
        print observed_tables
        assert expected_tables == observed_tables
    
    def count1(self,featuretype):
        self.c.execute('select count() from features where featuretype = ?',(featuretype,))
        results = self.c.fetchone()[0]
        print 'count1("%s") says: %s' % (featuretype,results)
        return results

    def count2(self,featuretype):
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
    
    def count3(self,featuretype):
        results = self.G.count_features_of_type(featuretype)
        print 'count3("%s") says: %s' % (featuretype, results)
        return results
    
    def count4(self,featuretype):
        cnt = 0
        for i in self.G.features_of_type(featuretype):
            cnt += 1
        print 'count4("%s") says: %s' % (featuretype,cnt)
        return cnt
    
    def featurecount_test(self):
        """check the number of each featuretype, using 4 different ways of counting."""
        featuretypes = ['gene',
                        'mRNA',
                        'CDS',
                        'exon',
                        'five_prime_UTR',
                        'three_prime_UTR',
                        'pcr_product',
                        'protein','intron']
        for featuretype in featuretypes:
            rawsql_cnt = self.count1(featuretype)

            # This is not an appropriate test for GTF files, since "gene" is
            # not explicitly listed -- only implied by the boundaries of CDSs:
            #
            # gffparsed_cnt = self.count2(featuretype)
            #
            count_feature_of_type_cnt = self.count3(featuretype)
            iterator_cnt = self.count4(featuretype)
            assert rawsql_cnt == count_feature_of_type_cnt == iterator_cnt

    def total_features_test(self):
        """did all features in the GFF file make it into the database?"""
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
        "parents test"
        if self.featureclass == 'GFF':
            parents1 = GFF_parent_check_level_1
            parents2 = GFF_parent_check_level_2
        if self.featureclass == 'GTF':
            parents1 = GTF_parent_check_level_1
            parents2 = GTF_parent_check_level_2

        for child, expected_parents in parents1.items():
            observed_parents = self.G.parents(child, level=1)
            for observed_parent in observed_parents:
                assert observed_parent.id in expected_parents

    def _parents(self):
        "parents test"

        # The two transcripts have the same gene parent.
        p1 = list(self.G.parents('FBtr0300690'))
        p2 = list(self.G.parents('FBtr0300689'))
        print p1
        print p2
        assert len(p1) == 1
        assert len(p2) == 1
        assert p1[0].id == p2[0].id

        # Make sure that the inverse is true.
        children_mRNAs = list(self.G.children(p1[0].id,featuretype='mRNA'))
        print children_mRNAs
        assert len(children_mRNAs) == 2 
        
        c = ['FBtr0300690', 'FBtr0300689']
        for i in children_mRNAs:
            assert i.id in c


        # now spot-check for some second-level parents.
        # Use a mix of exons, introns, UTRs, and CDSs.
        grandchildren = ['unnamed_exon_1', 'FBgn0031208:3', 'CDS_FBgn0031208:2_737', 'intron_FBgn0031208:2_FBgn0031208:4']
        
        for grandchild in grandchildren:
            grandparents = list(self.G.parents(grandchild,level=2))
            assert len(grandparents) == 1
            assert grandparents[0].id == 'FBgn0031208'

    def teardown(self):
        pass
        os.remove(testdbfn)

class TestGTFDBClass(TestGFFDBClass):
    featureclass = 'GTF'

