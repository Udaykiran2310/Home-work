# Home-work
# PubMed Papers Analyzer

A command-line tool to fetch research papers from PubMed and identify those with authors affiliated with pharmaceutical or biotech companies.

## Features

- Search PubMed using its full query syntax
- Identify papers with authors from pharmaceutical/biotech companies
- Export results to CSV or display in a formatted table
- Support for debug output and customizable result limits

## Installation

1. Make sure you have Python 3.8 or later installed
2. Install Poetry if you haven't already:
   ```bash
   pip install poetry
   ```
3. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/pubmed-papers-analyzer.git
   cd pubmed-papers-analyzer
   ```
4. Install dependencies:
   ```bash
   poetry install
   ```

## Usage

The tool can be used either through Poetry or after installation:

```bash
# Using Poetry
poetry run get-papers-list "your search query" [options]

# After installation
get-papers-list "your search query" [options]
```

### Options

- `query`: PubMed search query (required)
- `-f, --file`: Output file path (CSV format)
- `-d, --debug`: Enable debug output
- `-n, --max-results`: Maximum number of results to fetch (default: 100)
- `-e, --email`: Email address for PubMed API access (required)
- `-h, --help`: Show help message

### Example

```bash
poetry run get-papers-list "cancer therapy" -f results.csv -n 50 -e your.email@example.com
```

## Code Organization

The project is organized as follows:

```
pubmed-papers-analyzer/
├── src/
│   └── pubmed_papers/
│       ├── __init__.py
│       ├── core.py        # Core functionality and API interaction
│       └── cli.py         # Command-line interface
├── tests/                 # Test files
├── pyproject.toml        # Project configuration and dependencies
└── README.md            # This file
```

## Tools and Libraries Used

- [Poetry](https://python-poetry.org/): Dependency management and packaging
- [Biopython](https://biopython.org/): PubMed API interaction
- [Typer](https://typer.tiangolo.com/): Command-line interface
- [Rich](https://rich.readthedocs.io/): Terminal formatting and output
- [Pydantic](https://pydantic-docs.helpmanual.io/): Data validation and settings management

## Development

1. Install development dependencies:
   ```bash
   poetry install
   ```

2. Run tests:
   ```bash
   poetry run pytest
   ```

3. Format code:
   ```bash
   poetry run black src tests
   poetry run isort src tests
   ```

4. Type checking:
   ```bash
   poetry run mypy src
   ```

## License

This project is licensed under the MIT License. 
