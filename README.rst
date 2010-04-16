GFFutils
========
.. contents::

Overview
--------
This module is used for doing things with GFF files that are too
complicated for a simple ``awk`` or ``grep`` command line call.

For example, to get a BED file of genes from a GFF file, you can use something
simple like::

    grep "gene" chr2.gff | awk '{print $1,$4,$5}'

But how would you use commandline tools to get a BED file of 3' exons from
genes longer than 5 kb?  Or how would you get the average number of isoforms
for genes on the plus strand?  These more complex questions are actually quite
easy to answer using ``GFFutils``. 

A by-product of structuring a GFF file in this way is that it is easy to
convert to something like a refFlat format -- use the ``GFFDB.refFlat()``
method for this.

See the **Examples** below to jump right in, done, or follow along for a
step-by-step introduction.

Installation
------------

Unzip the source code, and from the source directory run (with root
priveliges)::
    
    python setup.py install

Now you're ready to create a GFF database and interact with it from a
Python shell like IPython.

Creating a GFFDB
----------------
For each GFF file you would like to use you need to create a GFF database.
This database is stored in a file, and is simply a sqlite3 database.  You
only have to do this once for each GFF file.  As long as you don't delete
the database from your hard drive, you don't have to do this time-consuming
step again.

The database will take roughly twice as much hard drive space as the
original text file.  This is the cost of interactivity and fast lookups.

You will need a GFF file to import.  From either a Python script or a Python
shell::

    import GFFutils
    
    # downloaded from, e.g., FlyBase
    gff_filename = '/data/annotations/dm3.gff'
    
    # the database that will be created
    db_filename = '/data/dm3.db'
    
    # do it! (this will take several minutes to run)
    GFFutils.create_gffdb(gfffn, db_filename)

Now ``dm3.db`` is the sqlite3 database that can be used in all sorts of
weird and wonderful ways, outlined below.

You can find more information on exactly what's going on in the "Strategy"
section below.

For power users, you can of course work on the database directly. Here's the
schema::

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


Using the database interactively
--------------------------------
First, wrap your new database in a ``GFFDB`` object::

    import GFFutils

    # Set up a GFFDB object, telling it the filename of your database as 
    # created above
    G = GFFutils.GFFDB('dm3.db')
    
.. note::
   
    For performance, most of the ``GFFDB`` class methods return iterators.  In
    practice, you will need to either convert them to a list or iterate through
    them in a list comprehension or a for-loop.  You can also grab the next item
    in an iterator with its ``.next()`` method.  All four ways of getting info
    from an iterator are shown below in the examples.

Identifying what's in the database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
What sorts of features are in the db?  The ``GFFDB.features`` method
returns an iterator of the featuretypes that were in the GFF file (and
which are now in the ``featuretype`` field of the sqlite3 database, which
this method accesses).  

Most methods in a ``GFFDB`` object return iterators for performance.

Since this is the first example of using the iterators returned by a
``GFFDB`` object, here are a few different ways to get the results
from the iterator it returns.
   
Method 0: Convert iterator to a list::

    featuretype_iterator = G.features()
    featuretypes = list(featuretype_iterator)

Method 1: Use iterator in a for-loop (preferred)::

    featuretype_iterator = G.features()
    for featuretype in featuretype_iterator:
        print featuretype

Method 2: Call ``next()`` incrementally on the iterator::

    featuretype_iterator = G.features()
    featuretype_1 = featuretype_iterator.next()
    featuretype_2 = featuretype_iterator.next()
    featuretype_3 = featuretype_iterator.next()
    featuretype_4 = featuretype_iterator.next()
    ...
    ...

    featuretypes = [featuretype1, featuretype2, ...]

It's mostly a matter of preference which method you use.  However, using
the for-loop approach is most memory-efficient, since only a single
featuretype is in memory at one time.  This is not too important for
iterating through featuretypes (of which there are usually <50; typically
3-10).  But when you want to iterate through 15,000 genes it can be useful.

In any case, we get something like the following.  What you see on your screen
depends entirely on the GFF file that you created your database from::
    
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
an iterator of ``GFFFeature`` objects.  These objects are described in
more detail in another section below.

``'gene'`` was in the list of ``featuretypes`` above.  Let's find out how many
genes there were. In this method, we're not bringing ALL the genes into a giant
list -- we'll just increment a counter.  Only a single ``GFFFeature`` object is
in memory at a time, which is the advantage of iterators . . . ::

    gene_count = 0
    for gene in G.features_of_type('gene'):
        gene_count += 1
    print gene_count
    
This is something I found myself doing quite often, so there's a shortcut method
that just does a ``count()`` in the SQL directly.  Use it like this::

    gene_count = G.count_features_of_type('gene')

Feature types not found in the db will not return an error (maybe
they should, eventually?); they just don't return anything::

    ncabbages = G.count_features_of_type('cabbage')
    print ncabbages  # zero cabbages.

Already know the ID of a feature?  Get the ``GFFFeature`` object
for that gene directly like this::

    my_favorite_feature = G['FBgn0002121']

I found myself getting a gene to play around with by doing this::

    g = G.features_of_type('gene').next()

However, this always returns the same gene.  For better testing, there's a
``random_feature()`` method that chooses a random feature out of the database.
You can specify a featuretype if you'd like; otherwise you have a chance of
getting any feature that was in the GFF file::

    g = G.random_feature('gene')

GFFFeatures in more detail
--------------------------
Just to make sure we're on the same page, here's the setup for this
section::

    import GFFutils
    G = GFFutils.GFFDB('dm3.db')

Let's get a single ``GFFFeature`` to work with::

    gene = G.random_feature('gene')

``GFFFeature`` objects, when printed, show useful information::

    GFFFeature gene 'FBgn0031208': chr2L:7529-9484 (+)
    #           ^          ^              ^         ^ 
    #           |          |              |         |
    # featuretype      accession   genomic coords   strand

``GFFFeature`` objects have an attribute, ``id``, which contains the
accession in the attributes field of the original GFF file::

    print gene.id

    'FBgn0031208'

They also have many other properties::

    print gene.start
    print gene.stop
    print gene.chr
    print gene.featuretype
    print gene.strand


You can get the length of a gene with::

    gene_len = gene.stop - gene.start

or you can use the perhaps-more-convenient::

    gene_len = len(gene)

In a ``GFFFeature`` object, the ``GFFFeature.attributes`` 
attribute holds all the info that was in the attributes column of your GFF
file.  This will vary based on what was in your original GFF file.  You can
get a list of this with::
    
    print gene.attributes._attrs

and you can access any of the attributes with a dot, then the attribute name.
For example, in the GFF file I used, the above code returned the following
available attributes::

    ['ID', 'Name', 'Ontology_term', 'Dbxref', 'derived_computed_cyto', 'gbunit']

So we could get the ontology terms for this gene with::

    print gene.attributes.Ontology_term

which shows something like::

    ['SO:0000010', 'SO:0000087', 'GO:0008234', 'GO:0006508']   

Or the DBxref (database cross-reference) for the gene with::

    print gene.attributes.Dbxref

  
You now know enough to be able to generate a line for a BED-format file (note
subtracting 1 from the start to convert to BED format's zero-based start)::

    line = '%s\t%s\t%s\t%s\t%s\t%s\n' % (gene.chr, 
                                         gene.start-1, 
                                         gene.stop, 
                                         gene.id, 
                                         gene.value, 
                                         gene.strand)
    print line

But ``GFFFeature`` objects have a convenience function,
``to_bed()``, which also accepts a number from 3 to 6 so you can tell it
how many BED fields you want returned (3 fields is the default).

So you could write a BED file of all the genes like so::

    fout = open('genes.bed','w')  # open a file for writing
    for i in G.features_of_type('gene'):
        fout.write(i.to_bed())
    fout.close()

This can be extremely useful for downstream processing by, for example,
BEDtools.


Navigating the hierarchy of features
------------------------------------

Here's how to find the transcripts belonging to a gene.  The
``GFFFeature.children`` and ``GFFFeature.parents()`` methods need a
feature ID as an argument, which is stored in the :attr:`GFFFeature.id`
attribute::

    for i in G.children(gene.id):
        print i

Here's how to find the exons belonging to a gene.  By default, level=1,
which means a 'hierarchy distance' of 1 (direct parent/children).  level=2
is analagous to grandparent/grandchild, which is used for the relationship
between genes/exons.  level=3 not currently implemented (not clear where it
would be used)::

    for i in G.children(gene_name, level=2):
        print i

Note that, depending on your GFF file, you may have more than just exons as
the children of genes (e.g., 3' UTRs, introns, 5' UTRs).  If you just want
the exons, then you can filter by feature type::

    for i in G.children(gene.id, level=2):
        if i.featuretype == 'exon':
            print i

File format conversions
-----------------------

Converting features to BED files was described above; briefly::

    fout = open('genes.bed','w')
    for gene in G.features_of_type('gene'):
        fout.write(gene.to_bed())
    fout.close()

Exporting a refFlat entry for one gene::

    print G.refFlat(gene_name)

Create a new file, writing a refFlat entry for each gene.  Note that the
``refFlat()`` method is set up such that it will return ``None`` if there
were no CDSs for a particular gene.  We don't want to write these to file,
but do want to keep track of them.

This will take a few seconds to run::
    
    missing_cds = []
    fout = open('mydatabase.refFlat','w')
    for gene in G.features_of_type('gene'):
        rflt = G.refFlat(gene.id)
        if rflt is not None:
            fout.write(rflt)
        else:
            missing_cds.append(gene)

    fout.close()

So, what were those genes that didn't have CDSs?  Check the first 25::
    
    for g in missing_cds[:25]:
        print g.attributes.Name[0]

Ahhhhh . . . a bunch of snoRNAs, tRNAs, etc.  Makes sense!


``GFFFeatures`` have a ``GFFFeature.tostring()`` method which prints
back the GFF file entry as a string (with the newline included).  This
makes it very easy to write new GFF files containing a subset of the
features in the original GFF file::

    # new GFF file with genes > 5kb
    fout = open('big-genes.gff','w')
    for gene in G.features_of_type('gene'):
        if len(gene) < 5000:
            fout.write(gene.tostring())
    fout.close()
    

Examples
--------

In each case, assume the following setup::

    import GFFutils
    GFFutils.create_gffdb('dm3.gff','dm3.db')
    G = GFFutils.GFFDB('dm3.db')

Gene count
~~~~~~~~~~

::

    gene_count = 0
    for gene in G.features_of_type('three_prime_UTR'):
        gene_count += 1
    print gene_count


Average gene length
~~~~~~~~~~~~~~~~~~~
::

    gene_lengths = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        gene_lengths += len(gene)
        gene_count += 1
    mean_gene_length = float(gene_lengths) / gene_count
    print mean_gene_length

Longest gene
~~~~~~~~~~~~
::

    maxlen = 0
    for gene in G.features_of_type('gene'):
        gene_len = len(gene)
        if gene_len > maxlen:
            maxlen = gene_len
    print maxlen

Average exon count
~~~~~~~~~~~~~~~~~~
::

    exon_count = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        gene_exon_count = 0
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
   lengths = [i.stop-i.start for i in G.features_of_type('exon')]
   p.hist(lengths,bins=50)
   p.show()


Average number of isoforms for genes on plus strand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    isoform_count = 0
    gene_count = 0
    for gene in G.features_of_type('gene'):
        if gene.strand == '-':
            continue
        isoforms = [i for i in G.children(gene.id) if i.featuretype=='mRNA']
        isoform_count += len(isoforms)
        gene_count += 1
    mean_isoform_count = float(isoform_count) / gene_count

            
Strategy
--------
A GFF database is built in several passes.  

During the first pass, the lines from the GFF file are split up into fields and 
imported into the ``features`` table.  If a "Parent" attribute is defined for the
feature, then we know its first-order parent and we can enter this into the ``relations`` 
table.

For example, say we have the following GFF line::

    chr2L FlyBase exon 8668 9276 .  + 0 ID=exon_1;Parent=mRNA_1

It will be entered into the ``features`` table like this

====== ===== ======= =========== ==== ===== ===== ====== ===== =================
ID     chrom source  featuretype start stop value strand phase attributes
====== ===== ======= =========== ==== ===== ===== ====== ===== =================
exon_1 chr2L FlyBase CDS         8668 9276  .     +      0     ID=exon_1;mRNA_1
====== ===== ======= =========== ==== ===== ===== ====== ===== =================

Since this CDS has an annotated parent, this relationship is entered into the ``relations`` table:

======= ======= =====
parent  child   level
======= ======= =====
mRNA_1  exon_1  1
======= ======= =====

Note that we can't assign any second-order parents.  On this first pass, we can
only add first-order parents because that's the only information that's
available on a single line in the GFF file.

At some point in the GFF file though, the parent transcript is found::

    chr2L FlyBase mRNA 7529 9484 . + . ID=mRNA_1;Parent=gene_1

...and we import it into the ``features`` table:

====== ===== ======= =========== ==== ===== ===== ====== ===== =======================
ID     chrom source  featuretype start stop value strand phase attributes
====== ===== ======= =========== ==== ===== ===== ====== ===== =======================
exon_1 chr2L FlyBase CDS         8668 9276  .     +      0     ID=exon_1;mRNA_1
mRNA_1 chr2L FlyBase mRNA        7529 9484  .     +      .     ID=mRNA_1;Parent=gene_1
====== ===== ======= =========== ==== ===== ===== ====== ===== =======================

as well as the ``relations`` table:

======= ======= =====
parent  child   level
======= ======= =====
mRNA_1  exon_1  1
gene_1  mRNA_1  1
======= ======= =====

...and these tables continue to grow as the GFF file is parsed.  When this
first pass is done, indexes are created to speed up searching in the second
pass.

The second pass looks at the ``relations`` table.  Note that **the current implementation
only goes 2 levels deep;** I still need to write a more general recursive form
of this to support hierarchies of arbitrary depth.

In the second pass, we go through each ID in the ``features`` column.  Then we
find that ID in the ``child`` column of the ``relations`` table and get its corresponding
parent.  In the example above, we find exon_1 in the ``child`` column.  Then we get its
parent (mRNA_1).  Then we take that parent and get *it's* parent by looking for it in the ``child`` column 
and then grabbing its parent (gene_1).

Now we know that gene_1 is the "grandparent" of exon_1, and we can enter it
into the ``relations`` table as a parent of level 2:

======= ======= =====
parent  child   level
======= ======= =====
mRNA_1  exon_1  1
gene_1  mRNA_1  1
gene_1  exon_1  2
======= ======= =====

In practice, the results of the "parent search" are written to a temporary text
file and then imported into the ``relations`` table as a batch in the end.
This is to avoid recalculating the index each time a new row is added, somthing
that would be extraordinarily time consuming.

Once the second pass is complete, indexes are built and the database is ready for use.

For a 130MB GFF file with 800,000+ features, the entire process takes a little
under 10 mins to run.  Luckily, you only need to make this time investment when
you have a new GFF file; if you already have a database built then using
GFFutils is quite fast.
