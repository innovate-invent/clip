from setuptools import setup, find_packages

with open('README.rst') as readme:
    setup(
        name='clipoverlap',
        version='1.0.8',
        packages=find_packages(),
        long_description=readme.read(),
        install_requires=['pysam', 'CigarIterator'],
        scripts=['bin/clip'],
        url='https://github.com/innovate-invent/clip',
        license='MIT',
        author='Nolan',
        author_email='innovate.invent@gmail.com',
        description='Clips overlapping regions in read mates of SAM/BAM files.',
        include_package_data=True
    )
