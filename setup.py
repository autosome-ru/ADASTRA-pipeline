from setuptools import setup, find_packages

setup(
    name='adastra',
    packages=find_packages(),
    version='1.0.0',
    entry_points={
        'console_scripts': [
            'adastra = scripts.start_point:main'
        ],
    },
    author="Sergey Abramov, Alexandr Boytsov",
    install_requires=[
        'numpy>=1.18.0',
        'pandas>=1.0.4',
        'matplotlib>=3.2.1',
        'seaborn>=0.10.1',
        'docopt>=0.6.2',
        'requests>=2.24.0',
        'statsmodels>=0.11.1',
        'scipy>=1.4.1',
        'babachi>=1.5.6'
    ],
    python_requires='>=3.6',
)
