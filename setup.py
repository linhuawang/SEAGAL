from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='seagal',
    version='2.5.2',
    description='Spatial Enrichment Analysis of Gene Association using L-index',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/linhuawang/SEAGAL',
    author='Linhua Wang',
    author_email='linhua.wang1994@gmail.com',
    license='MIT',
    zip_safe=False,
    packages=find_packages(),
    install_requires=[
          'scanpy[leiden]', 'numpy', 'anndata',
          'pandas', 'matplotlib', 'seaborn', 'scipy',
          'squidpy', 'scikit-learn', 'libpysal', 'statsmodels'
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
)
