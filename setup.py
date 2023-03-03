# -*- coding: utf-8 -*-
from distutils.core import setup
setup(
  name = 'pyfluidproperties',
  packages = ['pyfluidproperties'],
  version = '0.5',      
  license='agpl-3.0',
  description = 'Fluid properties for common fluids',
  author = 'Christoffer Rappmann',                   
  author_email = 'christoffer.rappmann@gmail.com',   
  url = 'https://github.com/ChristofferRa/pyfluidfroperties', 
  download_url = 'https://github.com/ChristofferRa/pyfluidfroperties/archive/v_01.tar.gz',
  keywords = ['Fluid', 'Properties', 'IAPWSIF97', 'Water', 'H2O', 'Water properties', 'Steam tables', 'Engineering thermodynamics'],
  install_requires=['numpy',
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Engineers',    
    'Topic :: Engineering :: Fluid mechanics',
    'License :: Affero General Public License v3.0',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
  ],
)