from setuptools import setup, find_packages

setup(
    name="satellite_visibility",
    version="1.0.0",
    packages=find_packages(),
    install_requires=["numpy>=1.20.0", "matplotlib>=3.0.0"],
    extras_require={
        "test": ["pytest>=6.0"]
    },
    author="Ryazanov Daniil",
    description="Библиотека для оценки видимости космических аппаратов",
    python_requires=">=3.6",
)