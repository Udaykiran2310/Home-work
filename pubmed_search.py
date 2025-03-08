"""Script to fetch and analyze PubMed papers."""

import sys
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
import re
import csv
import json
import xml.etree.ElementTree as ET
import time

import typer
from rich.console import Console
from rich.table import Table
import requests

app = typer.Typer()
console = Console()

class PubMedFetcher:
    """Class to handle PubMed API interactions."""
    
    def __init__(self, email: str):
        """Initialize the PubMed fetcher with user email."""
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.academic_keywords = {'university', 'college', 'institute', 'school',
                                'academia', 'laboratory', 'hospital', 'clinic',
                                'medical center', 'centre'}
        self.company_keywords = {'inc', 'corp', 'ltd', 'limited', 'llc',
                               'pharmaceutical', 'biotech', 'therapeutics',
                               'bioscience', 'laboratories'}
    
    def _is_company_affiliation(self, affiliation: str) -> bool:
        """Check if an affiliation is from a company."""
        if not affiliation:
            return False
        affiliation_lower = affiliation.lower()
        
        # Check for academic keywords
        if any(kw in affiliation_lower for kw in self.academic_keywords):
            return False
            
        # Check for company keywords
        return any(kw in affiliation_lower for kw in self.company_keywords)
    
    def _extract_email(self, text: str) -> Optional[str]:
        """Extract email address from text."""
        if not text:
            return None
        email_pattern = r'[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}'
        match = re.search(email_pattern, text)
        return match.group(0) if match else None
    
    def _parse_publication_date(self, pub_date_elem) -> datetime:
        """Parse PubMed publication date into datetime object."""
        year = 1900
        month = 1
        day = 1
        
        year_elem = pub_date_elem.find(".//Year")
        month_elem = pub_date_elem.find(".//Month")
        day_elem = pub_date_elem.find(".//Day")
        
        if year_elem is not None:
            year = int(year_elem.text)
        if month_elem is not None:
            try:
                month = int(month_elem.text)
            except ValueError:
                # Handle month names
                month_names = {
                    'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
                    'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
                }
                month = month_names.get(month_elem.text[:3], 1)
        if day_elem is not None:
            day = int(day_elem.text)
            
        return datetime(year, month, day)
    
    def fetch_papers(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Fetch papers from PubMed based on the query.
        
        Args:
            query: PubMed search query
            max_results: Maximum number of results to return
            
        Returns:
            List of paper dictionaries
        """
        # Search PubMed
        search_params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "usehistory": "y",
            "retmode": "xml",
            "tool": "PubMedPapersAnalyzer",
            "email": self.email
        }
        
        search_url = f"{self.base_url}/esearch.fcgi"
        search_response = requests.get(search_url, params=search_params)
        search_response.raise_for_status()
        
        search_root = ET.fromstring(search_response.text)
        id_list = [id_elem.text for id_elem in search_root.findall(".//Id")]
        
        if not id_list:
            return []
        
        papers = []
        
        # Fetch details for papers in batches
        batch_size = 50
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i + batch_size]
            
            # Fetch paper details
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(batch_ids),
                "retmode": "xml",
                "tool": "PubMedPapersAnalyzer",
                "email": self.email
            }
            
            fetch_url = f"{self.base_url}/efetch.fcgi"
            fetch_response = requests.get(fetch_url, params=fetch_params)
            fetch_response.raise_for_status()
            
            root = ET.fromstring(fetch_response.text)
            
            for article in root.findall(".//PubmedArticle"):
                # Extract basic information
                pubmed_id = article.find(".//PMID").text
                title_elem = article.find(".//ArticleTitle")
                title = title_elem.text if title_elem is not None else ""
                pub_date = self._parse_publication_date(article.find(".//PubDate"))
                
                # Process authors and affiliations
                authors = []
                company_affiliations = set()
                corresponding_email = None
                
                for author in article.findall(".//Author"):
                    last_name = author.find("LastName")
                    fore_name = author.find("ForeName")
                    name = f"{last_name.text if last_name is not None else ''} {fore_name.text if fore_name is not None else ''}".strip()
                    
                    affiliation_elem = author.find(".//Affiliation")
                    affiliation = affiliation_elem.text if affiliation_elem is not None else ""
                    email = self._extract_email(affiliation)
                    
                    is_non_academic = self._is_company_affiliation(affiliation)
                    if is_non_academic:
                        company_affiliations.add(affiliation)
                    
                    author_data = {
                        "name": name,
                        "affiliation": affiliation,
                        "email": email,
                        "is_corresponding": False,
                        "is_non_academic": is_non_academic
                    }
                    
                    if email and not corresponding_email:
                        corresponding_email = email
                        author_data["is_corresponding"] = True
                    
                    if is_non_academic:
                        authors.append(author_data)
                
                if authors:  # Only include papers with non-academic authors
                    papers.append({
                        "pubmed_id": pubmed_id,
                        "title": title,
                        "publication_date": pub_date,
                        "non_academic_authors": authors,
                        "company_affiliations": list(company_affiliations),
                        "corresponding_author_email": corresponding_email
                    })
            
            # Be nice to the API
            time.sleep(0.5)
        
        return papers

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
