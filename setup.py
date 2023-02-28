# -*- coding: utf-8 -*-
from distutils.core import setup
setup(
  name = 'pyfluidfroperties',         # How you named your package folder (MyLib)
  packages = ['pyfluidfroperties'],   # Chose the same as "name"
  version = '0.5',      # Start with a small number and increase it with every change you make
  license='agpl-3.0',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Fluid properties for common fluids',   # Give a short description about your library
  author = 'Christoffer Rappmann',                   # Type in your name
  author_email = 'christoffer.rappmann@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/ChristofferRa/pyfluidfroperties',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/ChristofferRa/pyfluidfroperties/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['Fluid', 'Properties', 'IAPWSIF97', 'Water', 'H2O'],   # Keywords that define your package best
  install_requires=['numpy',
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Engineers',      # Define that your audience are developers
    'Topic :: Engineering :: Fluid mechanics',
    'License :: Affero General Public License v3.0',   # Again, pick a license
    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
  ],
)