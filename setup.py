from setuptools import setup, find_namespace_packages


def parse_meta_yaml(file_path: str) -> list:
    """
    Parses a `meta.yaml` file to extract package dependencies listed under the 'run' section.

    This function reads a `meta.yaml` file line by line, identifies the 'run' section, and extracts
    the dependencies listed under it. Dependencies are returned as a list of strings.

    Args:
        file_path (str): The path to the `meta.yaml` file to be parsed.

    Returns:
        list: A list of dependency strings extracted from the 'run' section of the `meta.yaml` file.
              Excludes the 'python' dependency if present.

    Example:
        Given a `meta.yaml` file with the following content:
            requirements:
              run:
                - numpy
                - pandas
                - python >=3.8

        The function will return:
            ['numpy', 'pandas']

    Notes:
        - The function assumes that the 'run' section is followed by a list of dependencies,
          each prefixed with '- '.
        - The function stops parsing when it encounters an empty line after the 'run' section.
        - The 'python' dependency is explicitly excluded from the result.
    """
    requirements = []
    with open(file_path, 'r') as file:
        flag = False
        for line in file:
            line = line.strip()
            if line.startswith("run:"):
                flag = True
            elif flag and '- ' in line:
                m = line.split()
                if m[1] != "python":
                    requirements.append(line[2:].strip())
            elif flag and not line:
                break  # end with empty line
    return requirements



# solve meta.yaml
requirements = parse_meta_yaml('meta.yaml')


setup(
    name = "longbow",
    package_dir={"": "src"},
    packages=find_namespace_packages(where="src"),
    package_data={
        "longbow.model": ["*.csv"],
        "longbow.module": ["*.py"],
    },
    version = "2.3.0",
    description = "A Python program for nanopore sequencing basecalling configuration prediction",
    author = "Jun Mencius",
    author_email = "zjmeng22@m.fudan.edu.cn",
    url = "https://github.com/JMencius/longbow",
    keywords = ["longbow", "ont", "configuration"],
    python_requires = ">=3.7",
    install_requires = requirements,
    extras_require = {
        "dev": ["pytest"],
        },
    entry_points={
    "console_scripts": [
        "longbow = longbow.longbow:main",
        ],
    },
)

