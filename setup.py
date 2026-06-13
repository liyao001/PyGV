from setuptools import setup
import versioneer

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="GenomeViewer",
    python_requires=">=3.8",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Li Yao",
    author_email="regulatorygenome@gmail.com",
    description="Python Genome Viewer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["pygv", "pygv.configs", "pygv.errors", "pygv.scalers", "pygv.tracks", "pygv.tracks.logomaker"],
    install_requires=[
        "pydantic>=2,<3",
    ],
    url="https://github.com/liyao001/PyGV",
    license="GPL",
)
