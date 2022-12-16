from setuptools import setup, find_packages
setup(name='ABCD',
      version='0.1.0',
    #   url='https://github.com/liyy2/ABCD/',
    #   maintainer='Gerstein Lab',
    #   license='MIT',
    #   description='Deep learning models for wearables data',
    #   keywords=[
    #       'ABCD',
    #       'wearables',
    #       'consumer-grade wearables',
    #       'deep-learning',
    #   ],
      packages=find_packages(),
      install_requires=[
          'tsai=0.3.4',
          'sktime=0.14.1',
      ],)
    #   python_requires='>=3.7,<3.11')