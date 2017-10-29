import setuptools

setuptools.setup(
    name='peaks',
    version='0.1.1',
    url='https://github.com/rgiessmann/peaks',

    author='Robert Giessmann',
    author_email='r.giessmann@tu-berlin.de',

    description='Data analysis pipeline for DNA footprinting data.',
    long_description=open('README.md').read(),

    packages=setuptools.find_packages(exclude=['test', 'test.*']),

    platforms='any',

    include_package_data=True,
    package_data = {
        'scripts' : ['*']
    },

    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'scipy',
        'statsmodels'
    ],

    entry_points={
        'console_scripts': [
            'annotate=peaks.annotate:main',
        ],
    },

    classifiers = [],
)
