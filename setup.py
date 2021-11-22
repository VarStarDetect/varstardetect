from distutils.core import setup
setup(
  name = 'varstardetect',         # How you named your package folder (MyLib)
  packages = ['varstardetect'],   # Chose the same as "name"
  version = '1.1.0',      # Start with a small number and increase it with every change you make
  license='gpl-3.0',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'TESS Variable Star Light Curve Fitter',   # Give a short description about your library
  author = 'Nicolas Carrizosa Arias, Jorge Perez Gonzalez and Andres Cadenas Blanco',                   # Type in your name
  author_email = 'varstardetect@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/VarStarDetect/varstardetect',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/VarStarDetect/varstardetect/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['Star','Astronomy', 'Star Detection', 'Detection'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'validators',
          'beautifulsoup4',
      ],
  classifiers=[
    'Programming Language :: Python :: 3'
    'Framework :: Matplotlib'
    'Topic :: Scientific/Engineering :: Mathematics'
    'Topic :: Scientific/Engineering :: Astronomy'
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    'Development Status :: 4 - Beta'
    'Intended Audience :: Science/Research'
    'Operating System :: OS Independent'
  ],
)