try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Axion Analysis',
    'author': 'Ana Malagon',
    'url': 'https://github.com/atmalagon/python_axion_analysis',
    'download_url': 'https://github.com/atmalagon/python_axion_analysis.git',
    'author_email': 'atmalagon@gmail.com',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['core'],
    'scripts': [],
    'name': 'python_axion_analysis'
}

setup(**config)
