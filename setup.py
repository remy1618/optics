from setuptools import setup

setup(name='optics',
      version='0.0.7',
      description='Transfer matrix simulation of 1D optical multilayers',
      url='http://github.com/remy1618/optics',
      author='Remy Ko',
      author_email='remy.ko@mail.utoronto.ca',
      license='MIT',
      packages=['optics'],
      include_package_data=True,
      python_requires='>=2.7, <4',
      install_requires=[
          'numpy>=1.7',
          'matplotlib>=1.1',
          'scikit-image>=0.14',
          'scipy >=1.2'
      ])