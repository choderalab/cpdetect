from setuptools import setup

setup(name='cpdetect',
      description='Bayesian Change point detection',
      url='https://github.com/choderalab/cpdetect',
      author='Chaya D. Stern',
      packages=['cpdetect'],
      install_requires=[
            'numpy',
            'pandas',
            'scipy'
      ])