from setuptools import setup, find_packages

setup( 
        name='GFFutils',
        version='0.5',
        test_suite = 'test',
        py_modules = ['GFFutils'],
        scripts=['scripts/GFF_cleaner.py'],
        author='Ryan Dale',
        url='none',
        author_email='dalerr@niddk.nih.gov'

)

