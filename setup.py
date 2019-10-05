from setuptools import setup
from setuptools import find_packages


def setup_package():
    from peces import __version__
    REQUIRES = ['numpy', 'astropy>=1.0', 'pandas', 'matplotlib']
    META_DATA = dict(
        name='peces',
        version=__version__,
        description='Image reduction for CAFOS instrument in CAHA',
        author='Enrique Galceran',
        author_email='egalcera@ucm.es',
        packages=find_packages('.'),
        entry_points={
            'console_scripts': [
                'peces = peces.ULTRON:main',
                'peces-version = peces.version:main'
            ],
        },
        setup_requires=['peces'],
        install_requires=REQUIRES,
        zip_safe=False
    )

    setup(**META_DATA)


if __name__ == '__main__':
    setup_package()
