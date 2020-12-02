from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(
    name='pyHF',
    version='0.0.1',
    description='simple HF py code',
    long_description=readme,
    author='Amit',
    author_email='amit_@outlook.in',
    url='https://github.com/ipcamit/pyHF',
    packages=find_packages(exclude=('tests', 'docs'))
)
