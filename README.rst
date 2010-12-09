GFFutils
========
.. contents::

Overview and motivation
-----------------------
This module is used for doing things with GFF files that are too
complicated for a simple ``awk`` or ``grep`` command line call.

For example, to get a BED file of genes from a GFF file, you can use something
simple like::

    grep "gene" chr2.gff | awk '{print $1,$4,$5}' > out.bed

But how would you use commandline tools to get all the exons of a gene?  Or
exons that are found in all isoforms of a gene?  Or a BED file of 3' exons from
genes longer than 5 kb?  How would you get the average number of isoforms for
genes on the plus strand?  These more complex questions are actually quite easy
to answer using ``GFFutils``. 

A by-product of structuring a GFF file in this way is that it is easy to
convert to something like a refFlat format -- use the ``GFFDB.refFlat()``
method for this.

See the **Examples** below to jump right in, or follow along for a
step-by-step introduction.

Installation
------------

Unzip the source code, and from the source directory run (with root
priveliges)::
    
    python setup.py install

Now you're ready to create a GFF database and interact with it from a
Python shell like IPython.

New to Python?  Start here:

http://wiki.python.org/moin/BeginnersGuide

and here:

http://showmedo.com/videotutorials/python



Preparing to use GFFutils
-------------------------
For each GFF file you would like to use you need to create a GFF database.
This database is stored in a file, and is simply a sqlite3 database.  You
only have to do this once for each GFF file.  As long as you don't delete
the database from your hard drive, you don't have to do this time-consuming
step again.

The database will take roughly twice as much hard drive space as the
original text file.  This is the cost of interactivity and fast lookups.

You will need a GFF file to work with. If you don't already have one, here's one
you can use (*Drosophila melanogaster*; chromosome 2L only; 42 MB):

ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.29_FB2010/gff/dmel-2L-r5.29.gff.gz

If you're using a FlyBase GFF file, you might want to take a look at the
function ``GFFutils.clean_gff()`` in order to filter out features that may not
be of interest.  Assuming you have a cleaned-up GFF file, here's how to create
a GFF database from either a Python script or a Python shell.  To do this you
simply specify the path to the GFF file and the path to the new database you'd
like to create::

    >>> import GFFutils
    >>> GFFutils.create_gffdb('/data/annotations/dm3.gff', '/data/dm3.db')

This may take some time to run, but you only have to do this step once for
every GFF file you use.  Now ``dm3.db`` is the sqlite3 database that can be
used in all sorts of weird and wonderful ways, outlined below.

You can find more information on exactly what's going on in the `Strategy`_
section below.

For power users, you can of course work on the sqlite3 database directly.
Here's the schema::

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
    CREATE INDEX startstrand on features(start,strand);
    CREATE INDEX stops on features(stop);
    CREATE INDEX stopstrand on features(stop,strand);

Here's how to to connect to the new database from Python.  First, wrap your new
database in a ``GFFDB`` object::

    >>> G = GFFutils.GFFDB('dm3.db')

From now on we'll be accessing the database using this new object, ``G``, which
is a ``GFFutils.GFFDB`` object.

The next couple of sections will take the form of a tutorial. If you're itching
to get your hands dirty, all the methods should be documented so you can
explore the object interactively.  You might want to peek at the examples
below, too.


Identifying what's in the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For this section, I'm assuming that you've created a GFF database and have
connected to it as described above.  I'm also assuming you named the ``GFFDB``
object ``G``.

I'm also not assuming much Python knowledge.  If this sounds overly pedantic to
you, feel free to jump right to the `Examples`_!

As an introduction to using the database, let's start with answering a simple
question: "What sorts of features are in the GFF file?"  To do this, we'll use
the ``features()`` method of the ``GFFDB`` object.  The ``GFFDB.features()``
method returns a generator of the featuretypes that were in the GFF file (and
which are now in the ``featuretype`` field of the sqlite3 database, which this
method accesses).

Most methods in a ``GFFDB`` object return generators for performance.

.. note::
   
    For performance, most of the ``GFFDB`` class methods return generators.  In
    practice, you will need to either convert them to a list or iterate through
    them in a list comprehension or a for-loop.  You can also grab the next
    item in an iterator with its ``.next()`` method.  All four ways of getting
    info from a generator object are shown below in the examples.

Since this is the first example of using the generators returned by a ``GFFDB``
object, here are a few different ways to get the results from the generator.

Method 0: Convert iterator to a list.  This is the most memory-intensive::

    >>> featuretype_iterator = G.features()
    >>> featuretypes = list(featuretype_iterator)

Method 1: Use iterator in a for-loop (preferred)::

    >>> featuretype_iterator = G.features()
    >>> for featuretype in featuretype_iterator:
    ...     print featuretype

Method 2: Call ``next()`` incrementally on the iterator.  This is the most
awkward, but may sometimes be useful::

    >>> featuretype_iterator = G.features()
    >>> featuretype_1 = featuretype_iterator.next()
    >>> featuretype_2 = featuretype_iterator.next()
    >>> featuretype_3 = featuretype_iterator.next()
    >>> featuretype_4 = featuretype_iterator.next()
    ...
    ...

    >>> featuretypes = [featuretype1, featuretype2, ...]

It's mostly a matter of preference which method you use.  However, using
the for-loop approach is most memory-efficient, since only a single
featuretype is in memory at one time.  This is not too important for
iterating through featuretypes (of which there are usually <50; typically
3-10).  But when you want to iterate through 15,000 genes it can be useful.

In any case, we get something like the following.  What you see on your screen
depends entirely on the GFF file that you created your database from::
    
    >>> print featuretypes

    ['BAC_cloned_genomic_insert',
     'CDS',
     'DNA_motif',
     'breakpoint',
     'chromosome_arm',
     'chromosome_band',
     'complex_substitution',
     'deletion',
     'enhancer',
     'exon',
     'five_prime_UTR',
     'gene',
     'insertion_site',
     'intron',
     ...
     ...
      'tRNA',
     'tandem_repeat',
     'three_prime_UTR',
     'transposable_element',
     'transposable_element_insertion_site',
     'uncharacterized_change_in_nucleotide_sequence']



Retrieving specific feature types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To retrieve just genes, just exons, or any other feature type that was in
the GFF file, use the ``GFFDB.features_of_type()`` method.  This will return
an iterator of ``GFFFeature`` objects.  ``GFFFeature`` objects are described in
more detail in another section below.

``'gene'`` was in the list of ``featuretypes`` above.  Let's find out how many
genes there were. In this method, we're not bringing ALL the genes into a giant
list -- we'll just increment a counter.  Only a single ``GFFFeature`` object is
in memory at a time, which is the advantage of iterators . . . ::

    >>> gene_count = 0
    >>> for gene in G.features_of_type('gene'):
    ...     gene_count += 1
    >>> print gene_count
    
This is something I found myself doing quite often, so there's a shortcut method
that just does a ``count()`` in the SQL directly.  Use it like this::

    >>> gene_count = G.count_features_of_type('gene')

Feature types not found in the db will not return an error (maybe
they should, eventually?); they just don't return anything::

    >>> ncabbages = G.count_features_of_type('cabbage')
    >>> print ncabbages  # zero cabbages.

Already know the ID of a feature?  Get the ``GFFFeature`` object
for that gene directly like this::

    >>> my_favorite_feature = G['FBgn0002121']

The ID of a feature can be hard to remember.  The name of a gene is often much
easier to search by.  However, GFF files are not consistent in how they store
the name of a gene (for example, FlyBase GFF files have one name stored in the
Name attribute, while other names may be stored in the Alias attribute).
Nevertheless, there's a way to get named genes if the name is somewhere in the
attributes field::

    >>> candidates = G.attribute_search('Rm62')
    >>> assert len(candidates) == 1
    >>> my_favorite_gene = candidates[0]

This searches the attributes of all features of genes for the text 'Rm62'; the
search is case-insensitive.  Note that you get a list as a return value; that's
because there may be more than one gene with that text in the attributes; it's
up to you to figure out if the search returned the results you expected.

I found myself getting a gene to play around with by doing this::

    >>> g = G.features_of_type('gene').next()

However, this always returns the same gene.  For better testing, there's a
``random_feature()`` method that chooses a random feature out of the database.
You can specify a featuretype if you'd like; otherwise you have a chance of
getting any feature that was in the GFF file::

    >>> g = G.random_feature('gene')

GFFFeatures in more detail
--------------------------
This section discusses ``GFFFeatures`` which are the things you get back when
you query the database for a feature.

Just to make sure we're on the same page, here's the setup for this
section, assuming you've created a GFF database called ``dm3.db``::

    >>> import GFFutils
    >>> G = GFFutils.GFFDB('dm3.db')

Let's get a single ``GFFFeature`` to work with::

    >>> gene = G.random_feature('gene')

``GFFFeature`` objects, when printed, show useful information::

    GFFFeature gene 'FBgn0031208': chr2L:7529-9484 (+)
    #           ^          ^              ^         ^ 
    #           |          |              |         |
    # featuretype      accession   genomic coords   strand

``GFFFeature`` objects have an attribute, ``id``, which contains the
accession in the attributes field of the original GFF file::

    >>> print gene.id
    'FBgn0031208'

If there was no unique ID in the original GFF file, then the ID will be the
feature type plus the genomic coords (for example, "gene:chr2L:125-1340").  

``GFFFeature`` objects have many other properties::

    >>> gene.start
    7529

    >>> gene.stop
    9484

    >>> gene.chrom
    'chr2L'
    
    >>> gene.featuretype
    'gene'

    >>> gene.strand
    '+'

.. note::
   
   "name" is an alias to the "chrom" attribute if you're not working with
   features that can be mapped to a chromosome.  So in the code shown here, you
   can use ``gene.name`` instead of ``gene.chrom`` if it makes more semantic
   sense to do so for your application.

You can get the length of a feature with::

    >>> gene_len = gene.stop - gene.start

or you can use the perhaps-more-convenient::

    >>> gene_len = len(gene)

In a ``GFFFeature`` object, the ``GFFFeature.attributes`` 
attribute holds all the info that was in the attributes column of your GFF
file.  This will vary based on what was in your original GFF file.  You can
get a list of attribute names for a feature with::
    
    >>> print gene.attributes._attrs

and you can access any of the attributes with a dot, then the attribute name.
For example, in the GFF file I used, since the above code returned the following
available attributes::

    ['ID', 'Name', 'Ontology_term', 'Dbxref', 'derived_computed_cyto', 'gbunit']

then we could get the ontology terms for this gene with::

    >>> gene.attributes.Ontology_term
    ['SO:0000010', 'SO:0000087', 'GO:0008234', 'GO:0006508']   

Or the DBxref (database cross-reference) for the gene with::

    >>> gene.attributes.Dbxref
    ['FlyBase:FBan0011023',
     'FlyBase_Annotation_IDs:CG11023',
     'GB_protein:ACZ94128',
     'GB_protein:AAO41164',
     'GB:AI944728',
     'GB:AJ564667',
     'GB_protein:CAD92822',
     'GB:BF495604',
     'UniProt/TrEMBL:Q6KEV3',
     'UniProt/TrEMBL:Q86BM6',
     'INTERPRO:IPR003653',
     'EntrezGene:33155',
     'BIOGRID:59420',
     'FlyAtlas:CG11023-RA',
     'GenomeRNAi_gene:33155']

  
You now know enough to be able to generate a line for a BED-format file (note
subtracting 1 from the start to convert to BED format's zero-based start)::

    >>>line = '%s\t%s\t%s\t%s\t%s\t%s\n' % (gene.chrom, 
    ...                                     gene.start-1, 
    ...                                     gene.stop, 
    ...                                     gene.id, 
    ...                                     gene.value, 
    ...                                     gene.strand)
    >>> print line

But ``GFFFeature`` objects have a convenience function,
``to_bed()``, which also accepts a number from 3 to 6 so you can tell it
how many BED fields you want returned (3 fields is the default).

So you could write a BED file of all the genes longer than 5 kb like so::

    >>> fout = open('genes.bed','w')  # open a file for writing
    >>> for gene in G.features_of_type('gene'):
    ...     if len(gene) > 5000:
    ...         fout.write(gene.to_bed())
    >>> fout.close()

Other useful things in ``GFFFeature`` objects:

Reconstruct the GFF line for this feature, and automatically add a newline::

    >>> feature.tostring()

Get the transcription start site of the feature.  Note that all features have a
``TSS`` property, not just genes.  It is simply the feature's start position if
it's on the "+" strand or the feature's stop position if it's on the "-"
strand::

    >>> feature.TSS

Get the midpoint of the feature::

    >>> feature.midpoint

See the `Examples`_ below for more info on this.

Navigating the hierarchy of features
------------------------------------
Here's how to find the transcripts belonging to a gene.  The
``GFFDB.children()`` and ``GFFDB.parents()`` methods take either a feature ID
as an argument or a ``GFFFeature`` object.  The return value is a generator of
features that are children of the feature::

    >>> for i in G.children(gene.id):
    ...     print i

Here's how to find the exons belonging to a gene.  By default, level=1, which
means a 'hierarchy distance' of 1 (direct parent/children, for example genes
and transcripts).  level=2 is analagous to grandparent/grandchild, which is
used for the relationship between genes/exons.  level=3 is not currently
implemented::

    >>> for i in G.children(gene_name, level=2):
    ...     print i

Note that, depending on your GFF file, you may have more than just exons as
the children of genes (e.g., 3' UTRs, introns, 5' UTRs).  If you just want
the exons, then you can filter by feature type by specifying the
``featuretype`` keyword argument to ``children()``::

    >>> for exon in G.children(gene.id, level=2, featuretype='exon'):
    ...     print exon

Similarly, you can get the parents with ``GFFDB.parents()``.  Here's how to get
what gene an exon belongs to::

    >>> exon = G.random_feature('exon')
    >>> for grandparent_gene in G.parents(exon, level=2, featuretype='gene'):
    ...     print grandparent_gene

File format conversions
-----------------------

Converting features to BED files was described above; briefly::

    >>> fout = open('genes.bed','w')
    >>> for gene in G.features_of_type('gene'):
    ...     fout.write(gene.to_bed())
    >>> fout.close()

Exporting a refFlat entry for one gene::

    >>> print G.refFlat(gene_name)

Now create a new file, writing a refFlat entry for each gene.  Note that the
``refFlat()`` method is set up such that it will return ``None`` if there
were no CDSs for a particular gene.  We don't want to write these to file,
but do want to keep track of them.

This will take a few seconds to run::
    
    >>> missing_cds = []
    >>> fout = open('mydatabase.refFlat','w')
    >>> for gene in G.features_of_type('gene'):
    ...     rflt = G.refFlat(gene.id)
    ...     if rflt is not None:
    ...         fout.write(rflt)
    ...     else:
    ...         missing_cds.append(gene)
    >>> fout.close()

So, what were those genes that didn't have CDSs?  Check the first 25::
    
    >>> for g in missing_cds[:25]:
    ...     print g.attributes.Name[0]

A bunch of snoRNAs, tRNAs, etc.

``GFFFeatures`` have a ``GFFFeature.tostring()`` method which prints
back the GFF file entry as a string (with the newline included).  This
makes it very easy to write new GFF files containing a subset of the
features in the original GFF file::

    # new GFF file with genes > 5kb
    >>> fout = open('big-genes.gff','w')
    >>> for gene in G.features_of_type('gene'):
    ...     if len(gene) < 5000:
    ...         fout.write(gene.tostring())
    >>> fout.close()
    

Examples
--------

In each case, assume the following setup::

    import GFFutils
    GFFutils.create_gffdb('dm3.gff','dm3.db')
    G = GFFutils.GFFDB('dm3.db')

.. warning::

   These need examples need to be converted to doctests for thorough testing .
   . . email me if any of these don't work.

Inspecting the database
~~~~~~~~~~~~~~~~~~~~~~~
::
  
    print G.chromosomes()

    print G.strands()

    print list(G.features())


Gene count
~~~~~~~~~~
::

    G.count_features_of_type('gene')

Average gene length
~~~~~~~~~~~~~~~~~~~
::

    gene_lengths = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        gene_lengths += len(gene)
        gene_count += 1
    mean_gene_length = float(gene_lengths) / gene_count


BED file of genes
~~~~~~~~~~~~~~~~~
This will use the ID of the gene for the BED "name" field.

::


    fout = open('genes.bed','w')
    for gene in G.features_of_type('gene'):
        fout.write(gene.to_bed(6))
    fout.close()


Longest gene
~~~~~~~~~~~~
::

    maxlen = 0
    for gene in G.features_of_type('gene'):
        gene_len = len(gene)
        if gene_len > maxlen:
            maxlen = gene_len
            maxgene = gene
    print maxlen
    print maxgene

Longest gene, raw SQL
~~~~~~~~~~~~~~~~~~~~~
This version runs faster because it only ever looks at the start and stop
columns as opposed to the above version, which returns a full GFFFeature object
for each gene::

    c = G.conn.cursor()
    c.execute('''
        SELECT (stop-start) as LEN, * 
        FROM features
        WHERE featuretype="gene"
        ORDER BY LEN DESC
    ''')
    results = c.fetchone()
    maxlen = results[0]
       

Average exon count
~~~~~~~~~~~~~~~~~~
This takes several seconds to run, but as far as I know it's not something that
can be done easily using grep or awk::

    exon_count = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        gene_exon_count = 0

        # get all grandchildren, only counting the exons
        for child in G.children(gene.id,2):
            if child.featuretype == 'exon':
                gene_exon_count += 1

        exon_count += gene_exon_count
        gene_count += 1
    mean_exon_count = float(exon_count) / gene_count
    print mean_exon_count


Histogram of exon lengths
~~~~~~~~~~~~~~~~~~~~~~~~~
(Assumes you have matplotlib installed)

::

   from matplotlib import pyplot as p
   lengths = [len(i) for i in G.features_of_type('exon')]
   p.hist(lengths,bins=50)
   p.xlabel('Length of exon')
   p.ylabel('Frequency')
   p.show()


Average number of isoforms for genes on plus strand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This assumes that transcripts are labeled as "mRNA" instead of "transcript" or
something::

    isoform_count = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        if gene.strand == '-':
            continue
        isoforms = [i for i in G.children(gene.id) if i.featuretype=='mRNA']
        isoform_count += len(isoforms)
        gene_count += 1
    mean_isoform_count = float(isoform_count) / gene_count


BED file of all exonic bases on chr2L
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    exons = G.features_of_type('exon', chrom='chr2L')
    merged_exons = G.merge_features(exons,ignore_strand=True)
    fout = open('out.bed','w')
    for i in merged_exons:
        fout.write(i.to_bed())
    fout.close()

Closest features
~~~~~~~~~~~~~~~~

Get the closest gene (ignoring the gene you supply) and how far away it is::

    g = G.random_feature('gene')
    distance, closest_id = G.closest_feature(g.chr, 
                                             g.start,
                                             featuretype='gene',
                                             ignore=g.id)

Get the closest upstream exon that belongs to a different gene from the one you
supply::

    g = G.random_feature('gene')
    child_exons = G.children(g.id, level=2, featuretype='exon')
    ignore = [exon.id for exon in child_exons]
    distance, closest_exon = G.closest_feature(g.chr,
                                               g.start,
                                               featuretype='exon',
                                               ignore=ignore,
                                               strand=g.strand,
                                               direction='upstream')

Overlapping features
~~~~~~~~~~~~~~~~~~~~

Get the exons in the first MB of chr2L that are on the plus strand::

    exons_of_interest = G.overlapping_features(chrom='chr2L',
                                               start=1,
                                               stop=1e6,
                                               featuretype='exon',
                                               strand='+',
                                               completely_within=True)


Merging features
~~~~~~~~~~~~~~~~
This is useful if you want to get a "meta-exon" feature that is all exons
together.  For example, say you have a gene with two isoforms, and you want to
merge the exons together to get merged exons to indicate the presence of an
exon in *any* isoform.  Graphically::

                    
    isoform 1: [[[[[[[[[-----[[[[[[[[------------[[[
                 exon1         exon2             exon3
    isoform 2:     [[[[[[------------------------[[[[[[
                    exon4                         exon5

    merge    : [[[[[[[[[[----[[[[[[[[------------[[[[[[
                merged1         merged2           merged3

Code::

    g = G.random_feature('gene')
    exons = G.children(g.id, level=2, featuretype='exon')
    merged_exons = G.merge_features(exons)

    # If you want to create a new GFF file...
    fout = open('new.gff','w')
    for merged_exon in merged_exons:
        fout.write(merged_exon.tostring())
    fout.close()


Imputing introns
~~~~~~~~~~~~~~~~
Sometimes a GFF file doesn't explicitly include introns as features.  You can
construct them using the ``interfeatures()`` method.  This is a pretty
barebones method, so you'll have to add your own IDs and featuretypes after you
have the introns created.

::

    g = G.random_feature('gene')
    exons = G.children(g.id, level=2, featuretype='exon')
    introns = list(G.interfeatures(exons))
    
    for i,intron in enumerate(introns):
        intron.featuretype='intron'
        intron.add_attribute('ID', '%s_intron:%s' % (g.id,i))


Promoter regions
~~~~~~~~~~~~~~~~

Promoter regions, 1kb upstream and downstream of a gene's TSS::

    g = G.random_feature('gene')
    promoter = G.promoter(g.id)
    g.TSS - promoter.start
    promoter.stop - g.TSS

Promoter region defined as 2kb upstream::

    g = G.random_feature('gene')
    promoter = G.promoter(g.id, dist=2000, bidirectional=False)
    g.TSS - promoter.start
    promoter.stop - g.TSS

Coding genes
~~~~~~~~~~~~
Useful for excluding tRNAs, rRNAs, etc . . . this returns a generator of all
genes that have a CDS annotated as a child of level 2::

    we_make_proteins = G.coding_genes()

Isoform counts
~~~~~~~~~~~~~~
Useful for getting constitutive exons (i.e. exons found in all isoforms of a
gene)::

    g = G.random_feature('gene')
    n_gene_isos = G.n_gene_isoforms(g.id)
    for exon in G.children(g.id,level=2,featuretype='exon'): 
        if G.n_exon_isoforms(exon.id) == n_gene_isos:
            print exon.id, 'is found in all isoforms of', g.id

Arbitrary SQL commands
~~~~~~~~~~~~~~~~~~~~~~
Note that this places a lot of trust in the user to not mess up the database!

Things at the beginning of chromosomes::

    c = G.conn.cursor()
    results = c.execute("""
    SELECT id FROM features WHERE start BETWEEN 1 AND 100
    """)
    results = list(results)

Manually creating relationships::

    c.execute("""
    INSERT INTO relations VALUES ('fake_parent', 'fake_child', 100)
    """)

Manually removing relationships::

    c.execute("""
    DELETE FROM relations WHERE parent='fake_parent'
    """)


Strategy
--------
The following is my reasoning for the design of this package.  I'd be
interested to hear any thoughts on this or ways to improve it.

.. note::

   I tried a directed acyclic graph implementation, which would normally be
   useful for a hierachical data structure, but making it persistent meant
   unpickling it -- which took too long to start up and create.  Once it's
   created, the database approach seems to be the fastest.

A GFF database is built in several passes.  

During the first pass, the lines from the GFF file are split up into fields and 
imported into the ``features`` table.  If a "Parent" attribute is defined for the
feature, then we know its first-order parent and we can enter this into the ``relations`` 
table.

For example, say we have the following GFF line::

    chr2L FlyBase exon 8668 9276 .  + 0 ID=exon_1;Parent=mRNA_1

It will be entered into the ``features`` table like this::

    ID     chrom source  type start stop  value strand phase attributes
    ------ ----- ------- ---- ----- ----- ----- ------ ----- -----------------------
    exon_1 chr2L FlyBase CDS  8668  9276  .     +      0     ID=exon_1;mRNA_1

Since this CDS has an annotated parent, this relationship is entered into the ``relations`` table::

    parent  child   level
    ------- ------- -----
    mRNA_1  exon_1  1

Note that we can't assign any second-order parents.  On this first pass, we can
only add first-order parents because that's the only information that's
available on a single line in the GFF file.

At some point in the GFF file though, the parent transcript is found.  Here it is::

    chr2L FlyBase mRNA 7529 9484 . + . ID=mRNA_1;Parent=gene_1

...and we import it into the ``features`` table, just as the exon feature was added::

    ID     chrom source  type start stop  value strand phase attributes
    ------ ----- ------- ---- ----- ----- ----- ------ ----- -----------------------
    exon_1 chr2L FlyBase CDS  8668  9276  .     +      0     ID=exon_1;mRNA_1
    mRNA_1 chr2L FlyBase mRNA 7529  9484  .     +      .     ID=mRNA_1;Parent=gene_1

as well as the ``relations`` table, again just as the exon feature was added.
Note however that the mRNA_1 is now in the child column.  This will become
important later ::

    parent  child   level
    ------- ------- -----
    mRNA_1  exon_1  1
    gene_1  mRNA_1  1

The ``features`` table and the ``relations`` table continue to grow as the GFF
file is parsed.  Still, only first-order children/parents are added. When this
first pass is done, indexes are created to speed up searching in the second
pass.

The second pass looks at the ``relations`` table.  Note that **the current implementation
only goes 2 levels deep;** I still need to write a more general recursive form
of this to support hierarchies of arbitrary depth.

In the second pass, we go through each ID in the ``features`` column, matching
up IDs that are in the ``child`` column with the same ID in the ``parent``
column.  In the example above, we find "exon_1" in the ``child`` column.  Then
we get its parent ("mRNA_1").  Then we take that parent and get *it's* parent
by looking for "mRNA_1" in the ``child`` column and then grabbing its parent
("gene_1").

Now we know that gene_1 is the "grandparent" of exon_1, and we can enter it
into the ``relations`` table as a parent of level 2::

    parent  child   level
    ------- ------- -----
    mRNA_1  exon_1  1
    gene_1  mRNA_1  1
    gene_1  exon_1  2

In practice, the results of the "parent search" are written to a temporary text
file and then imported into the ``relations`` table as a batch in the end.
This is to avoid recalculating the index each time a new row is added, something
that would be extraordinarily time consuming.

Once the second pass is complete, indexes are built and the database is ready for use.

For a 130MB GFF file with 800,000+ features, the entire process takes a little
under 10 mins to run.  Luckily, you only need to make this time investment when
you have a new GFF file; if you already have a database built then using
GFFutils is quite fast.
