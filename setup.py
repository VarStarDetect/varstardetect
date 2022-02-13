from distutils.core import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='varstardetect',         # How you named your package folder (MyLib)
    packages=['varstardetect'],   # Chose the same as "name"
    version='0.2.1.3',      # Start with a small number and increase it with every change you make
    # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    license='gpl-3.0',
   
    # Give a short description about your library
    description="TESS Variable Star Light Curve Fitter",
    # Type in your name
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    author='Nicolas Carrizosa Arias, Jorge Perez Gonzalez and Andres Cadenas Blanco',
    author_email='varstardetect@gmail.com',      # Type in your E-Mail
    # Provide either the link to your github or to your website
    url='https://github.com/VarStarDetect/varstardetect',
    # I explain this later on
    download_url='https://github.com/VarStarDetect/varstardetect/archive/refs/tags/1.1.10.tar.gz',
    # Keywords that define your package best
    keywords=['Star', 'Astronomy', 'Star Detection', 'Detection'],
    install_requires=[

        'numpy',
        'matplotlib',
        'astropy',
        'pandas',
        'lightkurve',

    ],
    package_data={'varstardetect': ['*.txt','*.md','Targets/*.csv','*.in']},
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'Framework :: Matplotlib',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
    ],
    python_requires=">=3.8",
)
