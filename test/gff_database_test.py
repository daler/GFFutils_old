import tempfile, sqlite3, os, sys

exe_path = globals().get("__file__") or sys.argv[0]
exe_path = os.path.dirname(exe_path)
#sys.path.append(os.path.abspath(os.path.join(exe_path,"..","lib")))
#sys.path.append(os.path.abspath(os.path.join(exe_path,"..")))

import GFFutils

testdbfn = 'testing.db'

class TestGFFDBClass(object):
    def setup(self):
        if os.path.exists(testdbfn):
            os.remove(testdbfn)
        self.gfffn = os.path.join(exe_path,'test4.gff')
        self.dbfn = testdbfn
        GFFutils.create_gffdb(self.gfffn, self.dbfn)
        self.conn = sqlite3.connect(self.dbfn)
        self.c = self.conn.cursor()
        self.G = GFFutils.GFFDB(self.dbfn)
    
    def table_existence_test(self):
        """do the right tables exist?"""
        self.c.execute('select name from sqlite_master')
        tables = [i[0] for i in self.c]
        tables.sort()
        table_list = ['childindex',
                          'features',
                          'ids',
                          'parentindex',
                          'relations', 
                          u'sqlite_autoindex_features_1',
                          u'sqlite_autoindex_relations_1',
                          'starts',
                          'stops',
                          ]
        for i in tables:
            assert i in table_list

    
    def features_contents_test(self):
        """Contents of features OK?"""
        self.c.execute('select * from features')
        results = list(self.c)
        assert results == [(u'gn1', u'chrFAKE', 1,  170,  u'+',  u'gene',  u'.',  u'FlyBase',  u'.',  u'ID=gn1'),
                           (u'tr1', u'chrFAKE', 1,  170,  u'+',  u'mRNA',  u'.',  u'FlyBase',  u'.',  u'ID=tr1;Parent=gn1'),
                           (u'ex1', u'chrFAKE', 60, 93,   u'+',  u'exon',  u'.',  u'FlyBase',  u'.',  u'ID=ex1;Parent=tr1'),
                           (u'ex2', u'chrFAKE', 100,130,  u'+',  u'exon',  u'.',  u'FlyBase',  u'.',  u'ID=ex2;Parent=tr1')]
    
    def relations_contents_test(self):
        "Contents of relations OK?"
        self.c.execute('select * from relations')
        results = list(self.c)
        results.sort()
        for i in results:
            print '%s %s %s' % i
        assert results == [(u'gn1', u'ex1', 2),
                           (u'gn1', u'ex2', 2),
                           (u'gn1', u'tr1', 1),
                           (u'tr1', u'ex1', 1),
                           (u'tr1', u'ex2', 1)]

    def access_features_test(self):
        'test access of db' 
        genes = list(self.G.features_of_type('gene'))
        assert len(genes) == 1
        assert genes[0].__class__.__name__ == 'GFFFeature'
        #assert self.G[genes[0].id] == genes[0]
        transcripts = list(self.G.features_of_type('mRNA'))
        assert len(transcripts) == 1
        exons = list(self.G.features_of_type('exon'))
        assert len(exons) == 2
        
    def access_children_test(self):
        'test access of children'
        kids = list(self.G.children('gn1',level=1))
        print kids
        assert len(kids) == 1
        grandkids = list(self.G.children('gn1',level=2))
        print grandkids
        assert len(grandkids) == 2

    def access_parents_test(self):
        'test access of parents'
        parents = list(self.G.parents('tr1',level=1))
        print parents
        assert len(parents)==1

        grandparents = list(self.G.parents('ex1',level=2))
        print grandparents
        assert len(grandparents) == 1
    def access_all_test(self):
        'test access of all features'
        all = list(self.G.all())
        assert len(all) == 4

    def teardown(self):
        os.remove(testdbfn)
        

gtftestdbfn = 'testingGTF.db'
class TestGTFClass(object):
    def setup(self):
        if os.path.exists(testdbfn):
            os.remove(testdbfn)
        self.gfffn = os.path.join(exe_path,'test4.gtf')
        self.dbfn = testdbfn
        GFFutils.create_gtfdb(self.gfffn, self.dbfn)
        self.conn = sqlite3.connect(self.dbfn)
        self.c = self.conn.cursor()
        self.G = GFFutils.GTFDB(self.dbfn)
    
    def table_existence_test(self):
        """do the right tables exist?"""
        self.c.execute('select name from sqlite_master')
        tables = [i[0] for i in self.c]
        tables.sort()
        table_list = ['childindex',
                          'features',
                          'ids',
                          'parentindex',
                          'relations', 
                          u'sqlite_autoindex_features_1',
                          u'sqlite_autoindex_relations_1',
                          'starts',
                          'stops',
                          ]
        for i in tables:
            assert i in table_list
    
    def features_contents_test(self):
        """Contents of features OK?"""
        self.c.execute('select * from features order by id asc')
        results = list(self.c)
        for i in results:
            print i
        assert results == [(u'CDS:tr1:1', u'chrFAKE', 60, 93,   u'+',  u'CDS',  u'.',  u'FlyBase',  u'.',  u'gene_id gn1; transcript_id tr1; exon_number 1'),
                           (u'CDS:tr1:2', u'chrFAKE', 100,130,  u'+',  u'CDS',  u'.',  u'FlyBase',  u'.',  u'gene_id gn1; transcript_id tr1; exon_number 2'),
                           (u'gn1', u'chrFAKE', 60,  130,  u'+',  u'gene',  None, None, None, None),
                           (u'tr1', u'chrFAKE', 60,  130,  u'+',  u'mRNA',  None, None, None, None)]
    
    def relations_contents_test(self):
        "Contents of relations OK?"
        self.c.execute('select * from relations order by parent asc')
        results = list(self.c)
        results.sort()
        for i in results:
            print i
        assert results == [(u'gn1', u'CDS:tr1:1', 2),
                           (u'gn1', u'CDS:tr1:2', 2),
                           (u'gn1','tr1',1),
                           (u'tr1', u'CDS:tr1:1', 1),
                           (u'tr1', u'CDS:tr1:2', 1)]

    def access_features_test(self):
        'test access of db' 
        genes = list(self.G.features_of_type('gene'))
        assert len(genes) == 1
        assert genes[0].__class__.__name__ == 'GTFFeature'
        #assert self.G[genes[0].id] == genes[0]
        transcripts = list(self.G.features_of_type('mRNA'))
        assert len(transcripts) == 1
        exons = list(self.G.features_of_type('CDS'))
        assert len(exons) == 2
        
    def access_children_test(self):
        'test access of children'
        kids = list(self.G.children('gn1',level=1))
        print kids
        assert len(kids) == 1
        grandkids = list(self.G.children('gn1',level=2))
        print grandkids
        assert len(grandkids) == 2

    def access_parents_test(self):
        'test access of parents'
        parents = list(self.G.parents('tr1',level=1))
        print parents
        assert len(parents)==1

        grandparents = list(self.G.parents('CDS:tr1:1',level=2))
        print grandparents
        assert len(grandparents) == 1
    def access_all_test(self):
        'test access of all features'
        all = list(self.G.all())
        assert len(all) == 4

    def teardown(self):
        os.remove(testdbfn)
        


