import setuptools
setuptools.setup(
        name='optrotvib',
        author='Andrew M. James',
        description='tools for vibrational corrections to optical rotation',
        version='0.0.0',
        license='MIT',
        packages=setuptools.find_packages(),
        package_data={'optrotvib': ['schema/*.json']},
        entry_points = {
            "console_scripts": [
                'orv.engine': 'optrotvib.engine.cli:main'
                ]
            }
        install_requires = [
            'pyyaml',
            'py-cpuinfo',
            'psutil',
            'lxml',
            'numpy'
            ]
        zip_safe=False
        )
