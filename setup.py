#create a setup.py file to install the package
from setuptools import setup, find_packages


setup(
    name="synEvo",
    version="0.1",
    description="RNA benchmarking tools and utilities.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'einops',
        'forgi',
        "ViennaRNA",
        "biopython",
        "dask",
        "dask[dataframe]",
        "nltk",
        "regex",
        "Distance",
        "distlib",
        "gdown",
        "google-auth",
        "google-auth-oauthlib",
        "google-pasta",
        "GraKeL",
        "logomaker",
        "matplotlib",
        "numpy",
        "pandas",
        "pyaml",
        "pynisher",
        "pytest",
        "PyYAML",
        "rotary_embedding_torch",
        "scikit-bio",
        "seaborn",
        "tensorboard==2.12.3",
        "tensorboard-data-server==0.7.0",
        "tensorflow==2.12.0",
        "tensorflow-estimator==2.12.0",
        "tensorflow-io-gcs-filesystem==0.32.0",
        "torch",
        "torchvision",
        "umap",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.9",
        # Add other relevant classifiers.
    ],
# Optional
    author='Frederic Runge',
    author_email='runget@cs.uni-freiburg.de',
    license='MIT',
    keywords='Synthetic Evolution',
    url='https://github.com/Rungetf/RnaBench.git'
)