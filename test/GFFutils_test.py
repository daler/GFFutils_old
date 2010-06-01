import GFFutils
import os
import sqlite3

testdbfn = 'testing.db'

class TestGFFDBClass(object):
    def setup(self):
        if os.path.exists(testdbfn):
            os.remove(testdbfn)
        self.gfffn = 'FBgn0031208.gff'
        self.dbfn = testdbfn
        GFFutils.create_gffdb(self.gfffn, self.dbfn)
        self.conn = sqlite3.connect(self.dbfn)
        self.c = self.conn.cursor()
        self.G = GFFutils.GFFDB(self.dbfn)

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
        for line in open(self.gfffn):
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
        return self.G.count_features_of_type(featuretype)
    
    def count4(self,featuretype):
        cnt = 0
        for i in self.G.features_of_type(featuretype):
            cnt += 1
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
            gffparsed_cnt = self.count2(featuretype)
            count_feature_of_type_cnt = self.count3(featuretype)
            iterator_cnt = self.count4(featuretype)
            assert rawsql_cnt == gffparsed_cnt == count_feature_of_type_cnt == iterator_cnt

    def total_features_test(self):
        """did all features in the GFF file make it into the database?"""
        cnt = 0
        for feature in GFFutils.GFFFile(self.gfffn):
            cnt+=1
        self.c.execute('select count() from features')
        total = self.c.fetchone()[0]
        assert cnt == total

    def parents_test(self):
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
        #os.remove(testdbfn)
