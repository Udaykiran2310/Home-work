"""Command-line interface for the PubMed papers analyzer."""

import sys
from pathlib import Path
from typing import Optional
import csv
from datetime import datetime
import json

import typer
from rich.console import Console
from rich.table import Table

from .core import PubMedFetcher, Paper

app = typer.Typer()
console = Console()

def format_paper_for_csv(paper: dict) -> dict:
    """Format a paper object for CSV output."""
    return {
        "PubmedID": paper["pubmed_id"],
        "Title": paper["title"],
        "Publication Date": paper["publication_date"].strftime("%Y-%m-%d"),
        "Non-academic Author(s)": "; ".join(a["name"] for a in paper["non_academic_authors"]),
        "Company Affiliation(s)": "; ".join(paper["company_affiliations"]),
        "Corresponding Author Email": paper.get("corresponding_author_email", "")
    }

def save_to_csv(papers: list, output_file: Path):
    """Save papers to a CSV file."""
    if not papers:
        console.print("[yellow]No papers found matching the criteria.[/yellow]")
        return

    fieldnames = [
        "PubmedID",
        "Title",
        "Publication Date",
        "Non-academic Author(s)",
        "Company Affiliation(s)",
        "Corresponding Author Email"
    ]

    with output_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for paper in papers:
            writer.writerow(format_paper_for_csv(paper))

def print_table(papers: list):
    """Print papers in a formatted table."""
    if not papers:
        console.print("[yellow]No papers found matching the criteria.[/yellow]")
        return

    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("PubmedID", style="dim")
    table.add_column("Title")
    table.add_column("Publication Date")
    table.add_column("Non-academic Author(s)")
    table.add_column("Company Affiliation(s)")
    table.add_column("Corresponding Author Email")

    for paper in papers:
        table.add_row(
            paper["pubmed_id"],
            paper["title"],
            paper["publication_date"].strftime("%Y-%m-%d"),
            "\n".join(a["name"] for a in paper["non_academic_authors"]),
            "\n".join(paper["company_affiliations"]),
            paper.get("corresponding_author_email", "")
        )

    console.print(table)

@app.command()
def main(
    query: str = typer.Argument(..., help="PubMed search query"),
    output: Optional[Path] = typer.Option(
        None,
        "--file", "-f",
        help="Output file path (CSV format)"
    ),
    debug: bool = typer.Option(
        False,
        "--debug", "-d",
        help="Enable debug output"
    ),
    max_results: int = typer.Option(
        100,
        "--max-results", "-n",
        help="Maximum number of results to fetch"
    ),
    email: str = typer.Option(
        ...,
        "--email", "-e",
        help="Email address for PubMed API access",
        prompt="Please enter your email address for PubMed API access"
    )
):
    """
    Fetch research papers from PubMed and identify those with company-affiliated authors.
    
    The program searches PubMed using the provided query and identifies papers that have
    at least one author affiliated with a pharmaceutical or biotech company.
    """
    try:
        if debug:
            console.print(f"[blue]Debug: Query: {query}[/blue]")
            console.print(f"[blue]Debug: Max results: {max_results}[/blue]")

        fetcher = PubMedFetcher(email=email)
        
        with console.status("[bold green]Fetching papers from PubMed..."):
            papers = fetcher.fetch_papers(query, max_results=max_results)

        if output:
            save_to_csv(papers, output)
            console.print(f"[green]Results saved to {output}[/green]")
        else:
            print_table(papers)

    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        if debug:
            console.print_exception()
        sys.exit(1)

if __name__ == "__main__":
    app() 
